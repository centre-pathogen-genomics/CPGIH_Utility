#!/bin/bash

# USAGE: illumina_genomesqc.sh names inputdirectory outputdirectory

NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3

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
    
mkdir -p ${OUTPUTDIR}/KRAKEN/
mkdir -p ${OUTPUTDIR}/SPADES/

while IFS=$'\t' read -r i j k || [[ -n "$i" ]]
do

    echo 'Starting fastp processing of sample' ${i}
    echo 'Using reads in' ${j} ${k}

    fastp \
        --in1 ${j} \
        --in2 ${k} \
        --out1 FASTP/"$i"_R1_paired.fastq.gz \
        --out2 FASTP/"$i"_R2_paired.fastq.gz \
        --detect_adapter_for_pe \
        --length_required 50 \
        --thread 20 \
        --html FASTP/"$i"_fastp.html \
        --json FASTP/"$i"_fastp.json

    
    echo 'Starting Kraken2 classification of sample' ${i}
    echo 'Using reads in' ${j} ${k}

    kraken2 \
        --use-mpa-style \
        --use-names \
        --threads 20 \
        --paired \
        --output ${OUTPUTDIR}/KRAKEN/${i}_output.tsv \
        --report ${OUTPUTDIR}/KRAKEN/${i}_report.tsv \
        FASTP/"$i"_R1_paired.fastq.gz \
        FASTP/"$i"_R2_paired.fastq.gz

    rm -f ${OUTPUTDIR}/KRAKEN/${i}_output.tsv

    # pull out the 10 most abundant species from the report
    awk -F'\t' '$1 ~ /s__/ {gsub(/^ +| +$/, "", $0); print $0}' \
        ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | \
            sort -t$'\t' -k2,2nr | \
                head -n 10 > ${OUTPUTDIR}/KRAKEN/${i}_report_top10species.tsv

    # extract species counts from report - will use these after loop in summary output
    grep s__ ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | sed 's,.*s__,,' > ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv

    echo 'Starting Spades assembly of sample' ${i}
    echo 'Using reads in' ${j} ${k}
    
    spades.py \
        --isolate \
        -1 FASTP/"$i"_R1_paired.fastq.gz \
        -2 FASTP/"$i"_R2_paired.fastq.gz \
        -o ${OUTPUTDIR}/SPADES/${i}/ \
        -t 20

    mv ${OUTPUTDIR}/SPADES/${i}/contigs.fasta ${OUTPUTDIR}/SPADES/${i}_contigs.fa
    mv ${OUTPUTDIR}/SPADES/${i}/assembly_graph_with_scaffolds.gfa ${OUTPUTDIR}/SPADES/${i}_assembly_graph_with_scaffolds.gfa
    mv ${OUTPUTDIR}/SPADES/${i}/spades.log ${OUTPUTDIR}/SPADES/${i}_spades.log

    rm -rf ${OUTPUTDIR}/SPADES/${i}/

done < ${OUTPUTDIR}/.temp_manifest_filtered

# summarising kraken2 species results
echo -e "file\tspecies1\tspecies2\tspecies3" > ${OUTPUTDIR}/KRAKEN/top3species.tsv
# Loop through each report file
for file in ${OUTPUTDIR}/KRAKEN/*_report_species.tsv
do

    sample=$(basename "$file" _report_species.tsv)

    gawk -v sample="$sample" -F'\t' '
        {

            sum += $2
            data[NR] = $1
            counts[NR] = $2

        }
        END {
        output = sample
        if (sum == 0) {
        
            output = output "\tNA (0.00%)\tNA (0.00%)\tNA (0.00%)"
            
            } else {

                n = asorti(counts, idx, "@val_num_desc")

                for (e = 1; e <= 3; e++) {
                    if (e <= n) {
                        orig = idx[e]
                        species_name = data[orig]
                        percent = (counts[orig] / sum) * 100
                        output = output sprintf("\t%s (%.2f%%)", species_name, percent)
                    } else {
                        output = output "\tNA (0.00%)"
                    }
                }
            }
            print output
        }
    ' "$file" >> ${OUTPUTDIR}/KRAKEN/top3species.tsv

done

echo 'Computing assembly stats'
seqkit stats -abT ${OUTPUTDIR}/SPADES/*_contigs.fa | \
    cut -f 1,4,5,13 | \
    sed 's,_contigs.fa,,' | \
    sed 's,num_seqs,contigs, ; s,sum_len,assembly_length, ; s,N50,assembly_N50,' > ${OUTPUTDIR}/assembly_stats.tsv

csvtk join -t --left-join --na 0 -f file ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv | \
    cut -f 1,2,3,4,5,6,7,9,10,11,12,13,14 | \
    csvtk mutate2 -t  -n mean_coverage -e ' $sum_len / $assembly_length ' | \
    sed 's,+Inf,NA,'> ${OUTPUTDIR}/summary.tsv

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2 

# print information about empty reads sets
if [ "$SAMPLESREMOVED" -gt 0 ]
then

    echo ''
    echo 'The following samples were not analysed due to empty read sets:'
    cat ${OUTPUTDIR}/.emptysamples
    echo ''

else

    echo ''
    echo 'All sample read sets are non-empty, all were retained for analysis'
    echo ''

fi > emptysamples.info
