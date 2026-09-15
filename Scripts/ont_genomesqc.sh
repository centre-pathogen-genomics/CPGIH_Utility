#!/bin/bash

# USAGE: ont_genomesqc.sh names inputdirectory outputdirectory

NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3

# NCBI per-species genome length stats, used to sanity-check the LRGE genome size estimate
GSIZEDB="${CONDA_PREFIX}/stats_genomelength_ncbi.tsv"

# fail if errors are detected - only using during QC
set -e

# ensure names file exits
if [ ! -f ${NAMES} ]
then

    echo "Sample Names Input Does Not Exist. Mission Aborted."
    exit 1

fi

# ensure the NCBI genome size database exists - attempt to download it if missing
if [ ! -f "${GSIZEDB}" ]
then

    echo "NCBI genome size database not found at ${GSIZEDB} - attempting to download it"

    if curl -fsSL -o "${GSIZEDB}.tmp" "https://zenodo.org/records/21278902/files/stats_genomelength_ncbi.tsv"
    then

        mv "${GSIZEDB}.tmp" "${GSIZEDB}"
        echo "Downloaded NCBI genome size database to ${GSIZEDB}"

    else

        rm -f "${GSIZEDB}.tmp"
        echo "Failed to download NCBI genome size database. Mission Aborted."
        exit 1

    fi

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
mkdir -p ${OUTPUTDIR}/ASSEMBLIES/

# permanent record of the (possibly NCBI-corrected) genome size estimates for use in the summary
echo -e "file\tlrge_genome_size" > ${OUTPUTDIR}/lrge_gsize.tsv

# per-sample record of which assembler the Autocycler pipeline selected
echo -e "file\tassembler" > ${OUTPUTDIR}/assembler.tsv

# per-sample record of whether the genome size came from LRGE or the NCBI species median
echo -e "file\tgsize_source" > ${OUTPUTDIR}/gsize_source.tsv

while IFS=$'\t' read -r i j || [[ -n "$i" ]]
do

    echo 'Estimating genome size for sample' ${i} 'with LRGE'

    GSIZE=$(lrge ${j})
    GSIZE=$(printf '%.0f' "${GSIZE}")            # --genome_size requires an integer

    echo 'Starting Kraken2 classification of sample' ${i}
    echo 'Using reads in' ${j}

    kraken2 \
        --use-mpa-style \
        --use-names \
        --threads 16 \
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

    # sanity-check the LRGE estimate against the NCBI median genome size for the
    # most abundant classified species; fall back to the NCBI median if LRGE looks wrong
    TOPSPECIES=$(sort -t$'\t' -k2,2nr ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv | head -n 1 | cut -f1)

    if [ -n "${TOPSPECIES}" ]
    then

        DBMEDIAN=$(awk -F'\t' -v sp="s__${TOPSPECIES}" '$1==sp {print $6; exit}' "${GSIZEDB}")

    else

        DBMEDIAN=""

    fi

    if [ -n "${DBMEDIAN}" ]
    then

        DBMEDIAN=$(printf '%.0f' "${DBMEDIAN}")
        PCTDIFF=$(awk -v a="${GSIZE}" -v b="${DBMEDIAN}" 'BEGIN{d=(a-b)/b*100; if (d<0) d=-d; print d}')
        OUTOFRANGE=$(awk -v p="${PCTDIFF}" 'BEGIN{print (p>15)?"yes":"no"}')

        if [ "${OUTOFRANGE}" = "yes" ]
        then

            echo "LRGE genome size estimate for ${i} (${GSIZE} bp) is more than 15% from the NCBI median for ${TOPSPECIES} (${DBMEDIAN} bp) - using NCBI median instead"
            GSIZE=${DBMEDIAN}
            echo -e "${i}\tncbi_median" >> ${OUTPUTDIR}/gsize_source.tsv

        else

            echo -e "${i}\tlrge" >> ${OUTPUTDIR}/gsize_source.tsv

        fi

    else

        echo "No NCBI genome size database entry found for '${TOPSPECIES}' (sample ${i}) - keeping unverified LRGE estimate"
        echo -e "${i}\tlrge_unverified" >> ${OUTPUTDIR}/gsize_source.tsv

    fi

    echo -e "${i}\t${GSIZE}" >> ${OUTPUTDIR}/lrge_gsize.tsv

    echo 'Starting Autocycler (Flye fallback) assembly of sample' ${i}
    echo 'Using reads in' ${j}

    # NOTE: the output directory must not already exist - do not create it here
    autocycler_and_flye.py \
        --read-type ont_r10 \
        --genome_size ${GSIZE} \
        --threads 16 \
        --jobs 4 \
        --seed 42 \
        ${j} \
        ${OUTPUTDIR}/ASSEMBLIES/${i}

    if [ -f ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.fasta ]
    then

        cp ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.fasta ${OUTPUTDIR}/ASSEMBLIES/${i}_assembly.fasta

        ASSEMBLER=$(grep 'Final assembly' ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.log | sed 's,.* ,,' | tr '[:upper:]' '[:lower:]')
        [ -z "${ASSEMBLER}" ] && ASSEMBLER='NA'
        echo -e "${i}\t${ASSEMBLER}" >> ${OUTPUTDIR}/assembler.tsv

        # keep the plasmid summary alongside the assemblies (not merged into summary.tsv)
        if [ -f ${OUTPUTDIR}/ASSEMBLIES/${i}/plassembler_summary.tsv ]
        then
            cp ${OUTPUTDIR}/ASSEMBLIES/${i}/plassembler_summary.tsv \
                ${OUTPUTDIR}/ASSEMBLIES/${i}_plassembler_summary.tsv
        fi

    else

        echo 'WARNING: no assembly produced for sample' ${i}
        echo -e "${i}\tNA" >> ${OUTPUTDIR}/assembler.tsv

    fi

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
seqkit stats -abT ${OUTPUTDIR}/ASSEMBLIES/*_assembly.fasta | \
    cut -f 1,4,5,13 | \
    sed 's,_assembly.fasta,,' | \
    sed 's,num_seqs,contigs, ; s,sum_len,assembly_length, ; s,N50,assembly_N50,' > ${OUTPUTDIR}/assembly_stats.tsv

csvtk join -t --left-join --na 0 -f file ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv \
    ${OUTPUTDIR}/assembler.tsv \
    ${OUTPUTDIR}/gsize_source.tsv | \
    gawk -F'\t' -v OFS='\t' '
        # first file: build a sample -> LRGE genome size lookup
        NR==FNR {

            if (FNR > 1) gsize[$1] = $2
            next

        }

        # second file (joined summary from stdin): compute new columns
        FNR==1 {

            for (c=1; c<=NF; c++) col[$c] = c
            print $0, "mean_coverage", "predicted_genome_size", "lrge_qc", "coverage_qc", "assembly_qc", "contig_qc", "species_qc"
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

            # LRGE QC: flag implausible genome size estimates
            lrge_qc = (gs < 1800000 || gs > 6500000) ? "FLAG" : "PASS"

            # ASSEMBLY QC: fail if assembly length is more than 15% off the LRGE estimate
            if (gs > 0) {

                pct_diff = (assembly_length - gs) / gs * 100
                if (pct_diff < 0) pct_diff = -pct_diff
                assembly_qc = (pct_diff > 15) ? "FAIL" : "PASS"

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

            # COVERAGE QC
            if (mean_cov == "NA") {
                coverage_qc = "NA"
            } else if (mean_cov < 35) {
                coverage_qc = "FAIL"
            } else if (mean_cov < 40) {
                coverage_qc = "FLAG"
            } else {
                coverage_qc = "PASS"
            }

            print $0, mean_cov, gs, lrge_qc, coverage_qc, assembly_qc, contig_qc, species_qc

        }
    ' ${OUTPUTDIR}/lrge_gsize.tsv - > ${OUTPUTDIR}/summary.tsv

echo 'Formatting summary report'

# select/reorder the summary.tsv columns (dropping lrge_qc), give them more informative
# headers, and save as .xlsx
SUMMARYCOLS_OLD='file,reads,sum_len,min_len,avg_len,max_len,N50,contigs,assembly_length,assembly_N50,species1,species2,species3,assembler,predicted_genome_size,gsize_source,mean_coverage,coverage_qc,assembly_qc,contig_qc,species_qc'
SUMMARYCOLS_NEW='sequencing ID,number of reads,sum length,minimum read length,average read length,maximum read length,read N50,number of contigs,assembly length,assembly N50,kraken top species call,kraken species call 2,kraken species call 3,final assembly tool used,predicted genome size,predicted genome size tool used,mean coverage (sum length / predicted genome size),depth QC (mean coverage >40x = PASS; 35-40x = FLAG; <35x = FAIL),assembly QC (assembly length = predicted genome size +/- 15% PASS/FAIL),contig QC (>30 contigs = FAIL; 11-30 contigs = FLAG; <=10 contigs = PASS),species QC (kraken top species call >80% = PASS; top 3 calls same genus <80% = FAIL; else FLAG)'

csvtk cut -t -f "${SUMMARYCOLS_OLD}" ${OUTPUTDIR}/summary.tsv | \
    csvtk rename -t -f "${SUMMARYCOLS_OLD}" -n "${SUMMARYCOLS_NEW}" | \
    csvtk csv2xlsx -t -f -o ${OUTPUTDIR}/summary.xlsx

rm -f ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered ${OUTPUTDIR}/.temp_paths1 ${OUTPUTDIR}/.temp_paths2
rm -f ${OUTPUTDIR}/.temp_manifest.tsv ${OUTPUTDIR}/.temp_paths
rm -f ${OUTPUTDIR}/lrge_gsize.tsv ${OUTPUTDIR}/assembler.tsv ${OUTPUTDIR}/gsize_source.tsv

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
