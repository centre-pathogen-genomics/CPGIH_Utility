#!/bin/bash

# USAGE: illumina_genomesqc.sh names inputdirectory outputdirectory

set -e
NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3

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

    echo 'Creating output directory'${OUTPUTDIR}
    mkdir -p ${OUTPUTDIR}/

fi

# make manifest file
while read i
do

    ls ${INPUTDIR}/${i}*_R1_*.fastq.gz

done < ${NAMES} >${OUTPUTDIR}/.temp_paths1

while read i
do

    ls ${INPUTDIR}/${i}*_R2_*.fastq.gz

done < ${NAMES} >${OUTPUTDIR}/.temp_paths2

paste ${NAMES} ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 > ${OUTPUTDIR}/.temp_manifest
rm -f ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 

# ensure all specified input fastq files exist
FASTQERROR='false'
while read i j k 
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

mkdir -p ${OUTPUTDIR}/KRAKEN/
mkdir -p ${OUTPUTDIR}/SHOVILL/

while read i j k
do

    echo 'Starting Kraken2 classification of sample' ${i}
    kraken2 \
        --use-mpa-style \
        --use-names \
        --db /home/mdu/resources/kraken2/gtdb_r214/128gb/ \
        --confidence 0.1 \
        --threads 20 \
        --paired \
        --output ${OUTPUTDIR}/KRAKEN/${i}_output.tsv \
        --report ${OUTPUTDIR}/KRAKEN/${i}_report.tsv \
        ${j} \
        ${k}

    rm -f ${OUTPUTDIR}/KRAKEN/${i}_output.tsv

    # pull out the 10 most abundant species from the report
    awk -F'\t' '$1 ~ /s__/ {gsub(/^ +| +$/, "", $0); print $0}' \
        ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | \
            sort -t$'\t' -k2,2nr | \
                head -n 10 > ${OUTPUTDIR}/KRAKEN/${i}_report_top10species.tsv

    echo 'Starting Flye assembly of sample' ${i}
    shovill \
        --R1 ${j} \
        --R2 ${k} \
        --outdir ${OUTPUTDIR}/SHOVILL/${i}/ \
        --minlen 1000 \
        --cpus 20

    mv ${OUTPUTDIR}/SHOVILL/${i}/contigs.fa ${OUTPUTDIR}/SHOVILL/${i}_contigs.fa
    rm -rf ${OUTPUTDIR}/SHOVILL/${i}/

done < ${OUTPUTDIR}/.temp_manifest

rm -f ${OUTPUTDIR}/.temp_manifest

seqkit stats \
    -abT ${OUTPUTDIR}/SHOVILL/*_contigs.fa | \
        cut -f 1,4,5,13 | \
            csvtk pretty -t > ${OUTPUTDIR}/assembly_stats.tsv

