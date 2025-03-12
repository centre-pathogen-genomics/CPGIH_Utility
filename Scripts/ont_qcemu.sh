#!/bin/bash

# USAGE: ont_qcemu.sh names inputdirectory outputdirectory

set -e
NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3

MINLEN=1400
MAXLEN=1700

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


# ensure all specified input fastq files exist
FASTQERROR='false'
while IFS= read -r i || [[ -n "$i" ]]
do

    if [ ! -f ${INPUTDIR}/${i}.fastq.gz ]
	then

		echo 'File' ${i} 'does not exist'
		FASTQERROR='true'

	fi

done < ${NAMES}

# exit if fastq files don't exist
if [ ${FASTQERROR} = 'true' ]
then

    exit 1

fi

# make manifest file
while IFS= read -r i || [[ -n "$i" ]]
do

    ls ${INPUTDIR}/${i}.fastq.gz 

done < ${NAMES} > ${OUTPUTDIR}/.paths

paste -d $'\t' ${NAMES} ${OUTPUTDIR}/.paths > ${OUTPUTDIR}/.manifest.tsv
rm -f ${OUTPUTDIR}/.paths

# START PIPELINE

echo 'All specified inputs look good, starting pipeline'
echo 'Creating and outputting to' ${OUTPUTDIR}

mkdir -p ${OUTPUTDIR}/FILTERED_FASTQ/

while IFS=$'\t' read -r i j || [[ -n "$i" ]]
do
    
    echo 'Starting read filtering of sample' ${i}
	python /home/cwwalsh/Scripts/DAMG/ONT-16S/utils/filter_fastq.py \
	   --input_file ${j} \
	   --output_file ${OUTPUTDIR}/FILTERED_FASTQ/${i}.fastq.gz \
	   --min_length ${MINLEN} \
	   --max_length ${MAXLEN}

done < ${OUTPUTDIR}/.manifest.tsv

/home/cwwalsh/Software/seqkit stats\
	-abT ${OUTPUTDIR}/FILTERED_FASTQ/*.fastq.gz > ${OUTPUTDIR}/seqkit_stats_filtered.txt

while IFS=$'\t' read -r i j || [[ -n "$i" ]]
do 

    echo 'Starting Emu classification of sample' ${i}
	emu abundance \
		--type map-ont \
		--db /home/cwwalsh/Databases/Emu/ \
		--keep-counts \
		--output-dir ${OUTPUTDIR}/EmuResults/ \
		--output-basename ${i} \
		--threads 10 \
		${OUTPUTDIR}/FILTERED_FASTQ/${i}.fastq.gz

done < ${OUTPUTDIR}/.manifest.tsv

rm -f ${OUTPUTDIR}/EmuResults/*_rel-abundance-threshold-0.0001.tsv

emu combine-outputs --split-tables ${OUTPUTDIR}/EmuResults/ species
emu combine-outputs --split-tables --counts ${OUTPUTDIR}/EmuResults/ species

rm -f ${OUTPUTDIR}/.manifest.tsv


