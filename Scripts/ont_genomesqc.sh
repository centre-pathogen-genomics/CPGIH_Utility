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

done < ${NAMES} > ${OUTPUTDIR}/.temp_paths

paste ${NAMES} ${OUTPUTDIR}/.temp_paths > ${OUTPUTDIR}/.temp_manifest

# START PIPELINE

echo 'All specified inputs look good, starting pipeline'

echo 'Computing FASTQ read stats'
seqkit stats -abT --infile-list ${OUTPUTDIR}/.temp_paths | \
    cut -f 1,4,6,7,8,13 | \
    sed 's,.fastq.gz,,' | \
    sed 's,num_seqs,reads,' > ${OUTPUTDIR}/.read_stats

# identify empty read sets and remove from analysis loop
awk -F '\t' '$2 == 0' ${OUTPUTDIR}/.read_stats | cut -f 1 > ${OUTPUTDIR}/.emptysamples

if [ -s ${OUTPUTDIR}/.emptysamples ]; then
    awk -F '\t' 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
        ${OUTPUTDIR}/.emptysamples ${OUTPUTDIR}/.temp_manifest > ${OUTPUTDIR}/.temp_manifest_filtered
else
    cp ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered
fi

# remove empty read sets from read stats file
if [ -s ${OUTPUTDIR}/.emptysamples ]; then
    awk -F '\t' 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
        ${OUTPUTDIR}/.emptysamples ${OUTPUTDIR}/.read_stats > ${OUTPUTDIR}/read_stats.tsv
else
    cp ${OUTPUTDIR}/.read_stats ${OUTPUTDIR}/read_stats.tsv
fi

# print information about empty reads sets
echo 'Removing the following samples from QC due to empty read sets:'
cat ${OUTPUTDIR}/.emptysamples

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

    # extract species counts from report - will use these after loop in summary output
    grep s__ ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | sed 's,.*s__,,' > ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv
    
    echo 'Starting Flye assembly of sample' ${i}

    flye \
        --nano-hq ${j} \
        -o ${OUTPUTDIR}/FLYE/${i}/ \
        -t 20

    mv ${OUTPUTDIR}/FLYE/${i}/assembly.fasta ${OUTPUTDIR}/FLYE/${i}_assembly.fasta

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
seqkit stats -abT ${OUTPUTDIR}/FLYE/*_assembly.fasta | \
    cut -f 1,4,5,13 | \
    sed 's,_assmebly.fasta,,' | \
    sed 's,num_seqs,contigs, ; s,sum_len,assembly_length, ; s,N50,assembly_N50,' > ${OUTPUTDIR}/assembly_stats.tsv

paste ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv | \
    cut -f 1,2,3,4,5,6,8,9,10,11,12,13 > ${OUTPUTDIR}/summary.tsv

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 

grep 'Mean coverage' ${OUTPUTDIR}/FLYE/*/flye.log | \
    sed 's,.*FLYE/,, ; s,/flye.log:,, ; s,Mean coverage:\t,,' > ${OUTPUTDIR}/coverage_stats.tsv

rm -f ${OUTPUTDIR}/.temp_manifest.tsv ${OUTPUTDIR}/.temp_paths
