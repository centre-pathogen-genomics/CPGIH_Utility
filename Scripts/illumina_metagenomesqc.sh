#!/bin/bash

# USAGE: illumina_metagenomesqc.sh names inputdirectory outputdirectory

NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3
HOST=$4

# fail if errors are detected - only using during QC
set -e

# ensure names file exits
if [ ! -f ${NAMES} ]
then

    echo "Sample Names Input Does Not Exist. Mission Aborted."
    exit 1

fi

# ensure input directory exists
if [ ! -d ${INPUTDIR} ]
then

    echo "Input Directory Does Not Exist. Mission Aborted."
    exit 1

fi

# ensure output directory doesn't exist
# if it doesn't, create it
if [ -d ${OUTPUTDIR} ]
then

    echo "Output Directory Already Exists"
    exit 1

    else

    echo 'Creating output directory' ${OUTPUTDIR}
    mkdir -p ${OUTPUTDIR}/

fi

# ensure host name is set correctly
if [[ "$HOST" = 'human' ]] || [[ "$HOST" = 'mouse' ]]
then

    echo "Host Specified:" $HOST

    else

    echo "Valid host not specified, skipping decontamination step"

fi

# make manifest file
while IFS= read -r i || [[ -n "$i" ]]
do

    ls ${INPUTDIR}/${i}*_R1_*.fastq.gz

done < ${NAMES} > ${OUTPUTDIR}/.temp_paths1

while IFS= read -r i || [[ -n "$i" ]]
do

    ls ${INPUTDIR}/${i}*_R2_*.fastq.gz

done < ${NAMES} > ${OUTPUTDIR}/.temp_paths2

paste -d $'\t' ${NAMES} ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 > ${OUTPUTDIR}/.temp_manifest

# ensure all specified input fastq files exist
FASTQERROR='false'
while IFS=$'\t' read -r i j k  || [[ -n "$i" ]]
do
    
    if [ ! -f ${j} ]
	then

		echo 'File' ${j} 'does not exist'
		FASTQERROR='true'

	fi

	if [ ! -f ${k} ]
	then

		echo 'File' ${k} 'does not exist'
		FASTQERROR='true'

	fi

done < ${OUTPUTDIR}/.temp_manifest

# exit if fastq files don't exist
if [ ${FASTQERROR} = 'true' ]
then

    exit 1

fi

# START PIPELINE

echo 'All specified inputs look good, starting pipeline'

# removing error handling behaviour
set +e

echo 'Computing FASTQ read stats'
seqkit stats -abT --infile-list ${OUTPUTDIR}/.temp_paths1 | \
    cut -f 1,4,5,6,7,8,13 | \
    sed 's,_S.*.fastq.gz,,' | \
    sed 's,num_seqs,readpairs,' > ${OUTPUTDIR}/.read_stats

# identify empty read sets and remove from analysis loop
awk -F '\t' '$2 == 0' ${OUTPUTDIR}/.read_stats | cut -f 1 > ${OUTPUTDIR}/.emptysamples

if [ -s ${OUTPUTDIR}/.emptysamples ]
then

    awk -F '\t' 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
        ${OUTPUTDIR}/.emptysamples ${OUTPUTDIR}/.temp_manifest > ${OUTPUTDIR}/.temp_manifest_filtered

else

    cp ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered

fi

# remove empty read sets from read stats file
if [ -s ${OUTPUTDIR}/.emptysamples ]
then

    awk -F '\t' 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
        ${OUTPUTDIR}/.emptysamples ${OUTPUTDIR}/.read_stats > ${OUTPUTDIR}/read_stats.tsv

else

    cp ${OUTPUTDIR}/.read_stats ${OUTPUTDIR}/read_stats.tsv

fi

# print information about empty reads sets
SAMPLESREMOVED=$(wc -l < "${OUTPUTDIR}/.emptysamples")
if [ "$SAMPLESREMOVED" -gt 0 ]
then

    echo ''
    echo 'Removing the following samples from QC due to empty read sets:'
    cat ${OUTPUTDIR}/.emptysamples
    echo ''

else

    echo ''
    echo 'All sample read sets are non-empty, retaining all for analysis'
    echo ''

fi
    
mkdir -p ${OUTPUTDIR}/FASTP/
mkdir -p ${OUTPUTDIR}/KRAKEN/

while IFS=$'\t' read -r i j k || [[ -n "$i" ]]
do

    echo 'Starting fastp processing of sample' ${i}
    echo 'Using reads in' ${j} ${k}

    fastp \
        --in1 ${j} \
        --in2 ${k} \
        --out1 ${OUTPUTDIR}/FASTP/"$i"_R1_paired.fastq.gz \
        --out2 ${OUTPUTDIR}/FASTP/"$i"_R2_paired.fastq.gz \
        --detect_adapter_for_pe \
        --length_required 50 \
        --thread 20 \
        --html ${OUTPUTDIR}/FASTP/"$i"_fastp.html \
        --json ${OUTPUTDIR}/FASTP/"$i"_fastp.json
    
    # take host name from input 
    # human, mouse, none
    if [[ "$HOST" = 'human' ]] ; then
        
        echo "Removing $HOST Data From Sample $i"
        
        hostile clean \
            --fastq1 ${OUTPUTDIR}/FASTP/"$i"_R1_paired.fastq.gz \
            --fastq2 ${OUTPUTDIR}/FASTP/"$i"_R2_paired.fastq.gz \
            --aligner bowtie2 \
            --index /home/cwwalsh/Databases/Hostile/human-t2t-hla-argos985-mycob140 \
            --output "$OUTPUTDIR"/FASTP/ \
            --threads 20

    elif [[ "$HOST" = 'mouse' ]]; then

        hostile clean \
            --fastq1 ${OUTPUTDIR}/FASTP/"$i"_R1_paired.fastq.gz \
            --fastq2 ${OUTPUTDIR}/FASTP/"$i"_R2_paired.fastq.gz \
            --aligner bowtie2 \
            --index /home/cwwalsh/Databases/Hostile/mouse-mm39 \
            --output "$OUTPUTDIR"/FASTP/ \
            --threads 20
        
    else

        echo 'Skipping host removal'
        ln -s ${OUTPUTDIR}/FASTP/"$i"_R1_paired.fastq.gz ${OUTPUTDIR}/FASTP/"$i"_R1_paired.clean_1.fastq.gz
        ln -s ${OUTPUTDIR}/FASTP/"$i"_R2_paired.fastq.gz ${OUTPUTDIR}/FASTP/"$i"_R2_paired.clean_2.fastq.gz
        
    fi
    
    echo 'Starting Kraken2 classification of sample' ${i}

    kraken2 \
        --use-mpa-style \
        --use-names \
        --threads 20 \
        --paired \
        --output /dev/null \
        --report ${OUTPUTDIR}/KRAKEN/${i}_report.tsv \
        ${OUTPUTDIR}/FASTP/"$i"_R1_paired.clean_1.fastq.gz \
        ${OUTPUTDIR}/FASTP/"$i"_R2_paired.clean_2.fastq.gz

    # bracken

done < ${OUTPUTDIR}/.temp_manifest_filtered

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 
