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

echo 'Computing FASTQ read stats'
seqkit stats -abT --infile-list ${OUTPUTDIR}/.temp_paths1 | \
    cut -f 1,4,6,7,8,13 | \
    sed 's,_S.*.fastq.gz,,' | \
    sed 's,num_seqs,readpairs,' > ${OUTPUTDIR}/.read_stats

# identify empty read sets and remove from analysis loop
awk -F '\t' '$2 == 0' ${OUTPUTDIR}/.read_stats | cut -f 1 > ${OUTPUTDIR}/.emptysamples
[ -s ${OUTPUTDIR}/.emptysamples ] && awk -F '\t' 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
    ${OUTPUTDIR}/.emptysamples ${OUTPUTDIR}/.temp_manifest || cat ${OUTPUTDIR}/.temp_manifest > ${OUTPUTDIR}/.temp_manifest_filtered

# remove empty read sets from read stats file
[ -s ${OUTPUTDIR}/.emptysamples ] && awk -F '\t' 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
    ${OUTPUTDIR}/.emptysamples ${OUTPUTDIR}/.read_stats || cat ${OUTPUTDIR}/.read_stats > ${OUTPUTDIR}/read_stats.tsv

# print information about empty reads sets
echo 'Removing the following samples from QC due to empty read sets:'
cat ${OUTPUTDIR}/.emptysamples
    
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

    # extract species counts from report - will use these after loop in summary output
    grep s__ ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | sed 's,.*s__,,' > ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv

    echo 'Starting Spades assembly of sample' ${i}
    spades.py \
        --isolate \
        -1 ${j} \
        -2 ${k} \
        -o ${OUTPUTDIR}/SPADES/${i}/ \
        -t 20

    mv ${OUTPUTDIR}/SPADES/${i}/contigs.fasta ${OUTPUTDIR}/SPADES/${i}_contigs.fa
    mv ${OUTPUTDIR}/SPADES/${i}/assembly_graph_with_scaffolds.gfa ${OUTPUTDIR}/SPADES/${i}_assembly_graph_with_scaffolds.gfa
    mv ${OUTPUTDIR}/SPADES/${i}/spades.log ${OUTPUTDIR}/SPADES/${i}_spades.log

    rm -rf ${OUTPUTDIR}/SPADES/${i}/

done < ${OUTPUTDIR}/.temp_manifest_filtered

# summarising kraken2 species results
echo -e "species1\tspecies2\tspecies3" > ${OUTPUTDIR}/KRAKEN/top3species.tsv
# Loop through each report file
for file in ${OUTPUTDIR}/KRAKEN/*_report_species.tsv; do
    awk -F'\t' '
        {
            sum += $2
            data[NR] = $1
            counts[NR] = $2
        }
        END {
            if (sum == 0) {
                # If total sum is zero, output NA (0.00%) for all columns
                print "NA (0.00%)\tNA (0.00%)\tNA (0.00%)"
            } else {
                # Sort indexes based on counts
                n = asort(counts, sorted_counts, "@val_num_desc")

                # Create formatted output with up to 3 species
                for (e = 1; e <= 3; e++) {
                    if (e <= n) {
                        species_name = data[e]
                        percent = (counts[e] / sum) * 100
                        printf "%s (%.2f%%)", species_name, percent
                    } else {
                        printf "NA (0.00%%)"
                    }
                    if (e < 3) printf "\t"
                }
                print ""
            }
        }
    ' "$file" >> ${OUTPUTDIR}/KRAKEN/top3species.tsv
done

echo 'Computing assembly stats'
seqkit stats -abT ${OUTPUTDIR}/SPADES/*_contigs.fa | \
    cut -f 1,4,5,13 | \
    sed 's,_contigs.fa,,' | \
    sed 's,num_seqs,contigs, ; s,sum_len,assembly_length, ; s,N50,assembly_N50,' > ${OUTPUTDIR}/assembly_stats.tsv

paste ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv | \
    cut -f 1,2,3,4,5,6,8,9,10,11,12,13 > ${OUTPUTDIR}/summary.tsv

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 
