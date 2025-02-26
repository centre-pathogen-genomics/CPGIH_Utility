#!/bin/bash

# USAGE: illumina_genomesqc.sh names inputdirectory outputdirectory

set -e
NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3

# save directory that scripts are running from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# make manifest file
while read i
do

    ls ${INPUTDIR}/${i}*_R1_*.fastq.gz

done < ${NAMES} > ${OUTPUTDIR}/.temp_paths1

while read i
do

    ls ${INPUTDIR}/${i}*_R2_*.fastq.gz

done < ${NAMES} > ${OUTPUTDIR}/.temp_paths2

paste ${NAMES} ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 > ${OUTPUTDIR}/.temp_manifest

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
mkdir -p ${OUTPUTDIR}/SPADES/

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

    # extract species counts from report - will use these after looop in summary output
    grep s__ ${OUTPUTDIR}/KRAKEN/${i}_report.tsv > ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv

    echo 'Starting Spades assembly of sample' ${i}
    spades.py \
        --isolate \
        -1 ${j} \
        -2 ${k} \
        -o ${OUTPUTDIR}/SPADES/${i}/ \
        -t 20

    mv ${OUTPUTDIR}/SPADES/${i}/contigs.fasta ${OUTPUTDIR}/SPADES/${i}_contigs.fa

done < ${OUTPUTDIR}/.temp_manifest

echo 'Computing FASTQ read stats'
seqkit stats -abT --infile-list ${OUTPUTDIR}/.temp_paths1 | \
    cut -f 1,4 | \
    sed 's,_S.*.fastq.gz,,' | \
    sed 's,num_seqs,num_reads,' > ${OUTPUTDIR}/read_stats.tsv

echo 'Computing assembly stats'
seqkit stats -abT ${OUTPUTDIR}/SPADES/*_contigs.fa | \
    cut -f 1,4,5,13 | \
    sed 's,_contigs.fa,,' | \
    sed 's,num_seqs,num_contigs, ; s,sum_len,sum_len_contigs, ; s,N50,N50_contigs,' > ${OUTPUTDIR}/assembly_stats.tsv

Rscript ${SCRIPT_DIR}/krakenreport_top3species.R ${OUTPUTDIR}/KRAKEN/

paste ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv | \
    cut -f 1,2,4,5,6,8,9,10 > ${OUTPUTDIR}/summary.tsv

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 
