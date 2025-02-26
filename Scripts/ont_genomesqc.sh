#!/bin/bash

# USAGE: ont_genomesqc.sh names inputdirectory outputdirectory

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

    echo 'Creating output directory' ${OUTPUTDIR}
    mkdir -p ${OUTPUTDIR}/

fi

# ensure all specified input fastq files exist
FASTQERROR='false'
while read i 
do

    if [ ! -f ${INPUTDIR}/"$i".fastq.gz ]
	then

		echo 'File' "$j" 'does not exist'
		FASTQERROR='true'

	fi

done < ${NAMES}

# exit if fastq files don't exist
if [ ${FASTQERROR} = 'true' ]
then

    exit 1

fi

# make manifest file
while read i
do

    ls ${INPUTDIR}/${i}.fastq.gz 

done < ${NAMES} > ${OUTPUTDIR}/temp_paths

paste ${NAMES} ${OUTPUTDIR}/temp_paths > ${OUTPUTDIR}/temp_manifest.tsv

# START PIPELINE

echo 'All specified inputs look good, starting pipeline'

mkdir -p ${OUTPUTDIR}/KRAKEN/
mkdir -p ${OUTPUTDIR}/FLYE/

while read i j
do

    echo 'Starting Kraken2 classification of sample' ${i}

    kraken2 \
        --use-mpa-style \
        --use-names \
        --confidence 0.1 \
        --threads 20 \
        --output ${OUTPUTDIR}/KRAKEN/${i}_output.tsv \
        --report ${OUTPUTDIR}/KRAKEN/${i}_report.tsv \
        ${j}

    rm -f ${OUTPUTDIR}/KRAKEN/${i}_output.tsv

    # pull out the 10 most abundant species from the report
    awk -F'\t' '$1 ~ /s__/ {gsub(/^ +| +$/, "", $0); print $0}' \
        ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | \
            sort -t$'\t' -k2,2nr | \
                head -n 10 > ${OUTPUTDIR}/KRAKEN/${i}_report_top10species.tsv

    echo 'Starting Flye assembly of sample' ${i}

    flye \
        --nano-hq ${j} \
        -o ${OUTPUTDIR}/FLYE/${i}/ \
        -t 20

    mv ${OUTPUTDIR}/FLYE/${i}/assembly.fasta ${OUTPUTDIR}/FLYE/${i}_assembly.fasta

done < ${OUTPUTDIR}/temp_manifest.tsv

/home/cwwalsh/Software/seqkit stats \
    -abT ${OUTPUTDIR}/FLYE/*_assembly.fasta | \
        cut -f 1,4,5,13 > ${OUTPUTDIR}/assembly_stats.tsv

grep 'Mean coverage' ${OUTPUTDIR}/FLYE/*/flye.log | \
    sed 's,.*FLYE/,, ; s,/flye.log:,, ; s,Mean coverage:\t,,' > ${OUTPUTDIR}/coverage_stats.tsv

rm -f ${OUTPUTDIR}/temp_manifest.tsv ${OUTPUTDIR}/temp_paths
