#!/bin/bash

# USAGE: ont_genomesqc.sh names inputdirectory outputdirectory

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

# ensure all specified input fastq files exist
FASTQERROR='false'
while IFS=$'\t' read -r i || [[ -n "$i" ]]
do

    if [ ! -f ${INPUTDIR}/"$i".fastq.gz ]
	then

		echo 'File' "$i" 'does not exist'
		FASTQERROR='true'

	fi

done < ${NAMES}

# exit if fastq files don't exist
if [ ${FASTQERROR} = 'true' ]
then

    exit 1

fi

# make manifest file
while IFS=$'\t' read -r i || [[ -n "$i" ]]
do

    ls ${INPUTDIR}/${i}.fastq.gz 

done < ${NAMES} > ${OUTPUTDIR}/.temp_paths

paste -d $'\t' ${NAMES} ${OUTPUTDIR}/.temp_paths > ${OUTPUTDIR}/.temp_manifest

# START PIPELINE

echo 'All specified inputs look good, starting pipeline'

# removing error handling behaviour
set +e

echo 'Computing FASTQ read stats'
seqkit stats -abT --infile-list ${OUTPUTDIR}/.temp_paths | \
    cut -f 1,4,5,6,7,8,13 | \
    sed 's,.fastq.gz,,' | \
    sed 's,num_seqs,reads,' > ${OUTPUTDIR}/.read_stats

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
mkdir -p ${OUTPUTDIR}/FLYE/

# subsample reads to 100x depth ahead of classification/assembly
mkdir -p ${OUTPUTDIR}/SUBSAMPLE/

while IFS=$'\t' read -r i j || [[ -n "$i" ]]
do

    echo 'Subsampling sample' ${i} 'to 100x depth'

    GSIZE=$(lrge ${j})
    echo -e "${i}\t${GSIZE}" >> ${OUTPUTDIR}/.temp_gsize
    rasusa reads -g ${GSIZE} -c 100 -s 42 ${j} -o ${OUTPUTDIR}/SUBSAMPLE/${i}_100x.fastq.gz

done < ${OUTPUTDIR}/.temp_manifest_filtered

# build a permanent record of the LRGE genome size estimates for use in the summary
echo -e "file\tlrge_genome_size" > ${OUTPUTDIR}/lrge_gsize.tsv
cat ${OUTPUTDIR}/.temp_gsize >> ${OUTPUTDIR}/lrge_gsize.tsv

# build manifest pointing to the subsampled reads
awk -F '\t' -v dir="${OUTPUTDIR}/SUBSAMPLE" '{print $1"\t"dir"/"$1"_100x.fastq.gz"}' \
    ${OUTPUTDIR}/.temp_manifest_filtered > ${OUTPUTDIR}/.temp_manifest_subsampled

while IFS=$'\t' read -r i j || [[ -n "$i" ]]
do

    echo 'Starting Kraken2 classification of sample' ${i}
    echo 'Using reads in' ${j}

    kraken2 \
        --use-mpa-style \
        --use-names \
        --threads 20 \
        --output /dev/null \
        --report ${OUTPUTDIR}/KRAKEN/${i}_report.tsv \
        ${j}

    # pull out the 10 most abundant species from the report
    awk -F'\t' '$1 ~ /s__/ {gsub(/^ +| +$/, "", $0); print $0}' \
        ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | \
            sort -t$'\t' -k2,2nr | \
                head -n 10 > ${OUTPUTDIR}/KRAKEN/${i}_report_top10species.tsv

    # extract species counts from report - will use these after loop in summary output
    grep s__ ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | sed 's,.*s__,,' > ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv
    
    echo 'Starting Flye assembly of sample' ${i}
    echo 'Using reads in' ${j}

    flye \
        --nano-hq ${j} \
        -o ${OUTPUTDIR}/FLYE/${i}/ \
        -t 20

    mv ${OUTPUTDIR}/FLYE/${i}/assembly.fasta ${OUTPUTDIR}/FLYE/${i}_assembly.fasta

done < ${OUTPUTDIR}/.temp_manifest_subsampled

# remove subsampled reads now that classification/assembly is complete
rm -rf ${OUTPUTDIR}/SUBSAMPLE/

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
seqkit stats -abT ${OUTPUTDIR}/FLYE/*_assembly.fasta | \
    cut -f 1,4,5,13 | \
    sed 's,_assembly.fasta,,' | \
    sed 's,num_seqs,contigs, ; s,sum_len,assembly_length, ; s,N50,assembly_N50,' > ${OUTPUTDIR}/assembly_stats.tsv

csvtk join -t --left-join --na 0 -f file ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv | \
    gawk -F'\t' -v OFS='\t' '
        # first file: build a sample -> LRGE genome size lookup
        NR==FNR {

            if (FNR > 1) gsize[$1] = $2
            next

        }

        # second file (joined summary from stdin): compute new columns
        FNR==1 {

            for (c=1; c<=NF; c++) col[$c] = c
            print $0, "mean_coverage", "lrge_genome_size", "assembler", "lrge_qc", "assembly_qc", "contig_qc", "species_qc"
            next

        }

        function species_name(str) {
            match(str, /^(.*) \([0-9.]+%\)$/, m)
            return m[1]
        }

        function species_pct(str) {
            match(str, /^.* \(([0-9.]+)%\)$/, m)
            return m[1] + 0
        }

        function species_genus(str,    parts) {
            split(str, parts, " ")
            return parts[1]
        }

        {

            sum_len = $(col["sum_len"])
            contigs = $(col["contigs"])
            assembly_length = $(col["assembly_length"])
            sp1 = $(col["species1"]); sp2 = $(col["species2"]); sp3 = $(col["species3"])

            gs = ($1 in gsize) ? gsize[$1] + 0 : 0

            # mean coverage: total read bases / LRGE predicted genome size
            mean_cov = (gs > 0) ? sum_len / gs : "NA"

            assembler = "flye"

            # LRGE QC: flag implausible genome size estimates
            lrge_qc = (gs < 1800000 || gs > 6500000) ? "FLAG" : "PASS"

            # ASSEMBLY QC: fail if assembly length is more than 10% off the LRGE estimate
            if (gs > 0) {

                pct_diff = (assembly_length - gs) / gs * 100
                if (pct_diff < 0) pct_diff = -pct_diff
                assembly_qc = (pct_diff > 10) ? "FAIL" : "PASS"

            } else {

                assembly_qc = "NA"

            }

            # CONTIG QC
            if (contigs > 30) {
                contig_qc = "FAIL"
            } else if (contigs > 10) {
                contig_qc = "FLAG"
            } else {
                contig_qc = "PASS"
            }

            # SPECIES QC: based on top3 species and the genus of the top hit
            name1 = species_name(sp1)

            if (name1 == "NA") {

                species_qc = "FAIL"

            } else {

                pct1 = species_pct(sp1)
                name2 = species_name(sp2); pct2 = species_pct(sp2)
                name3 = species_name(sp3); pct3 = species_pct(sp3)

                genus1 = species_genus(name1)
                genus_sum = pct1
                if (name2 != "NA" && species_genus(name2) == genus1) genus_sum += pct2
                if (name3 != "NA" && species_genus(name3) == genus1) genus_sum += pct3

                if (genus_sum < 80) {
                    species_qc = "FAIL"
                } else if (pct1 < 80) {
                    species_qc = "FLAG"
                } else {
                    species_qc = "PASS"
                }

            }

            print $0, mean_cov, gs, assembler, lrge_qc, assembly_qc, contig_qc, species_qc

        }
    ' ${OUTPUTDIR}/lrge_gsize.tsv - > ${OUTPUTDIR}/summary.tsv

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered ${OUTPUTDIR}/.temp_manifest_subsampled ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2
rm -f ${OUTPUTDIR}/.temp_manifest.tsv ${OUTPUTDIR}/.temp_paths ${OUTPUTDIR}/.temp_gsize
rm -f ${OUTPUTDIR}/lrge_gsize.tsv

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
