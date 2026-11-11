#!/bin/bash

# USAGE: ont_genomesqc.sh names inputdirectory outputdirectory

NAMES=$1
INPUTDIR=$2
OUTPUTDIR=$3

# NCBI per-species genome length stats, used to sanity-check the LRGE genome size estimate
# built by Scripts/build_ncbi_genome_size_db.py and stored in the repo's data/ directory
# set GSIZEDB in the environment to use a different file
SCRIPTDIR=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)
GSIZEDB="${GSIZEDB:-${SCRIPTDIR}/../data/prok_species_genome_size_stats.tsv}"

# fail if errors are detected - only using during QC
set -e

# ensure names file exits
if [ ! -f ${NAMES} ]
then

    echo "Sample Names Input Does Not Exist. Mission Aborted."
    exit 1

fi

# ensure the NCBI genome size database exists
if [ ! -f "${GSIZEDB}" ]
then

    echo "NCBI genome size database not found at ${GSIZEDB}. Mission Aborted."
    echo "Expected in the repo's data/ directory (or set GSIZEDB to its location)."
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
csvtk filter -t -f 'reads=0' ${OUTPUTDIR}/.read_stats | \
    csvtk cut -t -f file | \
    csvtk del-header -t > ${OUTPUTDIR}/.emptysamples

if [ -s ${OUTPUTDIR}/.emptysamples ]
then
    
    csvtk grep -t -H -f 1 -v -P ${OUTPUTDIR}/.emptysamples \
        ${OUTPUTDIR}/.temp_manifest > ${OUTPUTDIR}/.temp_manifest_filtered

else
   
   cp ${OUTPUTDIR}/.temp_manifest ${OUTPUTDIR}/.temp_manifest_filtered

fi

# remove empty read sets from read stats file
csvtk filter -t -f 'reads>0' ${OUTPUTDIR}/.read_stats > ${OUTPUTDIR}/read_stats.tsv

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
mkdir -p ${OUTPUTDIR}/LRGE/

# raw LRGE genome size estimates (before any NCBI correction) with 95% confidence intervals
echo -e "file\tlrge_estimate\tlrge_ci95_lower\tlrge_ci95_upper" > ${OUTPUTDIR}/lrge_estimates.tsv

# permanent record of the (possibly NCBI-corrected) genome size estimates for use in the summary
echo -e "file\tlrge_genome_size" > ${OUTPUTDIR}/lrge_gsize.tsv

# per-sample record of which assembler the Autocycler pipeline selected
echo -e "file\tassembler" > ${OUTPUTDIR}/assembler.tsv

# per-sample record of whether the genome size came from LRGE or the NCBI species median
echo -e "file\tgsize_source" > ${OUTPUTDIR}/gsize_source.tsv

while IFS=$'\t' read -r i j || [[ -n "$i" ]]
do

    echo 'Estimating genome size for sample' ${i} 'with LRGE'

    # estimate goes to stdout; the log (incl. the 95% CI) goes to stderr - keep it and echo it
    GSIZE=$(lrge ${j} 2> ${OUTPUTDIR}/LRGE/${i}_lrge.log)
    cat ${OUTPUTDIR}/LRGE/${i}_lrge.log >&2
    GSIZE=$(printf '%.0f' "${GSIZE}")            # --genome_size requires an integer

    # record the raw LRGE estimate and its 95% CI (in bp) before any NCBI correction
    # LRGE only reports the CI in its log, e.g. "(95% CI: 4.21 Mbp - 4.87 Mbp)", rounded to 2 dp
    CI=$(sed -n 's/.*95% CI: \([0-9.]*\) \([kMG]*bp\) - \([0-9.]*\) \([kMG]*bp\)).*/\1 \2 \3 \4/p' \
        ${OUTPUTDIR}/LRGE/${i}_lrge.log | head -n 1 | tr ' ' '\t')
    if [ -n "${CI}" ]
    then
        CI=$(echo -e "lo\tlo_unit\thi\thi_unit\n${CI}" | \
            csvtk mutate2 -t -n lower -w 0 \
                -e '$lo * ($lo_unit == "Gbp" ? 1000000000 : ($lo_unit == "Mbp" ? 1000000 : ($lo_unit == "kbp" ? 1000 : 1)))' | \
            csvtk mutate2 -t -n upper -w 0 \
                -e '$hi * ($hi_unit == "Gbp" ? 1000000000 : ($hi_unit == "Mbp" ? 1000000 : ($hi_unit == "kbp" ? 1000 : 1)))' | \
            csvtk cut -t -f lower,upper | \
            csvtk del-header -t)
    else
        CI=$'NA\tNA'
    fi
    echo -e "${i}\t${GSIZE}\t${CI}" >> ${OUTPUTDIR}/lrge_estimates.tsv

    echo 'Starting Kraken2 classification of sample' ${i}
    echo 'Using reads in' ${j}

    kraken2 \
        --use-mpa-style \
        --use-names \
        --threads 8 \
        --output /dev/null \
        --report ${OUTPUTDIR}/KRAKEN/${i}_report.tsv \
        ${j}

    # pull out the 10 most abundant species from the report
    csvtk grep -t -H -f 1 -r -p 's__' ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | \
        csvtk sort -t -H -k 2:nr | \
            head -n 10 > ${OUTPUTDIR}/KRAKEN/${i}_report_top10species.tsv

    # extract species counts from report - will use these after loop in summary output
    grep s__ ${OUTPUTDIR}/KRAKEN/${i}_report.tsv | sed 's,.*s__,,' > ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv

    # sanity-check the LRGE estimate against the NCBI median genome size for the
    # most abundant classified species; fall back to the NCBI median if LRGE looks wrong
    TOPSPECIES=$(sort -t$'\t' -k2,2nr ${OUTPUTDIR}/KRAKEN/${i}_report_species.tsv | head -n 1 | cut -f1)

    if [ -n "${TOPSPECIES}" ]
    then

        # prefer median of complete genomes; fall back to median of all genomes,
        # which is more affected by incomplete assemblies/MAGs (columns looked up by header name)
        DBMEDIAN=""
        DBMEDIANTYPE=""
        read -r NCOMPLETE MEDCOMPLETE MEDALL < <(
            csvtk grep -t -f species -p "${TOPSPECIES}" "${GSIZEDB}" | \
                csvtk cut -t -f n_complete,median_complete,median_all | \
                    csvtk del-header -t | \
                        head -n 1) || true   # read returns 1 if species not in DB

        if [ -n "${NCOMPLETE}" ] && [ "${NCOMPLETE}" != "NA" ] && [ "${NCOMPLETE}" -gt 0 ] && [ "${MEDCOMPLETE}" != "NA" ]
        then
            DBMEDIAN=${MEDCOMPLETE}; DBMEDIANTYPE='complete'
        elif [ -n "${MEDALL}" ] && [ "${MEDALL}" != "NA" ]
        then
            DBMEDIAN=${MEDALL}; DBMEDIANTYPE='all'
        fi

    else

        DBMEDIAN=""
        DBMEDIANTYPE=""

    fi

    if [ -n "${DBMEDIAN}" ]
    then

        DBMEDIAN=$(printf '%.0f' "${DBMEDIAN}")
        # integer arithmetic: |GSIZE - DBMEDIAN| / DBMEDIAN > 15%
        ABSDIFF=$(( GSIZE - DBMEDIAN )); ABSDIFF=${ABSDIFF#-}

        if (( ABSDIFF * 100 > DBMEDIAN * 15 ))
        then

            echo "LRGE genome size estimate for ${i} (${GSIZE} bp) is more than 15% from the NCBI median (${DBMEDIANTYPE} genomes) for ${TOPSPECIES} (${DBMEDIAN} bp) - using NCBI median instead"
            GSIZE=${DBMEDIAN}
            echo -e "${i}\tncbi_median_${DBMEDIANTYPE}" >> ${OUTPUTDIR}/gsize_source.tsv

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
        --read_type ont_r10 \
        --genome_size ${GSIZE} \
        --threads 8 \
        --jobs 4 \
        --seed 42 \
        ${j} \
        ${OUTPUTDIR}/ASSEMBLIES/${i}

    if [ -f ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.fasta ]
    then

        cp ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.fasta ${OUTPUTDIR}/ASSEMBLIES/${i}_assembly.fasta
        cp ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.gfa ${OUTPUTDIR}/ASSEMBLIES/${i}_assembly.gfa

        ASSEMBLER=$(grep 'Final assembly' ${OUTPUTDIR}/ASSEMBLIES/${i}/assembly.log | sed 's,.* ,,' | tr '[:upper:]' '[:lower:]')
        [ -z "${ASSEMBLER}" ] && ASSEMBLER='NA'
        echo -e "${i}\t${ASSEMBLER}" >> ${OUTPUTDIR}/assembler.tsv

        # keep the plasmid summary alongside the assemblies (not merged into summary.tsv)
        if [ -f ${OUTPUTDIR}/ASSEMBLIES/${i}/plassembler_summary.tsv ]
        then
            cp ${OUTPUTDIR}/ASSEMBLIES/${i}/plassembler_summary.tsv ${OUTPUTDIR}/ASSEMBLIES/${i}_plassembler_summary.tsv
        fi

    else

        echo 'WARNING: no assembly produced for sample' ${i}
        echo -e "${i}\tNA" >> ${OUTPUTDIR}/assembler.tsv

    fi

done < ${OUTPUTDIR}/.temp_manifest_filtered

# summarising kraken2 species results: top 3 species per sample as "species (percent%)",
# padded with "NA (0.00%)" when fewer than 3 species (or no reads) were classified
echo -e "file\tspecies1\tspecies2\tspecies3" > ${OUTPUTDIR}/KRAKEN/top3species.tsv
for file in ${OUTPUTDIR}/KRAKEN/*_report_species.tsv
do

    sample=$(basename "$file" _report_species.tsv)
    TOTAL=$(csvtk summary -t -H -f 2:sum "$file" 2>/dev/null || true)

    {
        if [ -n "${TOTAL}" ] && [ "${TOTAL}" != "0.00" ]
        then
            csvtk sort -t -H -k 2:nr "$file" | \
                head -n 3 | \
                    csvtk mutate2 -t -H -w 2 -e "\$2 / ${TOTAL} * 100" | \
                        csvtk replace -t -H -f 3 -p '^(.*)$' -r '($1%)' | \
                            csvtk cut -t -H -f 1,3 | \
                                tr '\t' ' '
        fi
        printf 'NA (0.00%%)\nNA (0.00%%)\nNA (0.00%%)\n'
    } | head -n 3 | paste -s -d $'\t' - | sed "s,^,${sample}\t," >> ${OUTPUTDIR}/KRAKEN/top3species.tsv

done

echo 'Computing assembly stats'
seqkit stats -abT ${OUTPUTDIR}/ASSEMBLIES/*_assembly.fasta | \
    cut -f 1,4,5,13 | \
    sed 's,_assembly.fasta,,' | \
    sed 's,num_seqs,contigs, ; s,sum_len,assembly_length, ; s,N50,assembly_N50,' > ${OUTPUTDIR}/assembly_stats.tsv

# join all per-sample tables, then compute mean coverage and the QC calls
#   lrge_qc:     FLAG if predicted genome size is outside 1.8-6.5 Mb
#   coverage_qc: mean coverage <35x FAIL; 35-40x FLAG; >=40x PASS
#   assembly_qc: FAIL if assembly length is more than 15% off the predicted genome size
#   contig_qc:   0 or >30 contigs FAIL; 11-30 FLAG; <=10 PASS
#   species_qc:  top-hit genus (summed over top 3 calls) <80% FAIL; top species <80% FLAG; else PASS
csvtk join -t --left-join --na 0 -f file ${OUTPUTDIR}/read_stats.tsv \
    ${OUTPUTDIR}/assembly_stats.tsv \
    ${OUTPUTDIR}/KRAKEN/top3species.tsv \
    ${OUTPUTDIR}/assembler.tsv \
    ${OUTPUTDIR}/gsize_source.tsv \
    ${OUTPUTDIR}/lrge_gsize.tsv | \
    csvtk rename -t -f lrge_genome_size -n predicted_genome_size | \
    csvtk mutate2 -t -n mean_coverage -w 2 \
        -e '$predicted_genome_size > 0 ? $sum_len / $predicted_genome_size : "NA"' | \
    csvtk mutate2 -t -n lrge_qc \
        -e '$predicted_genome_size < 1800000 || $predicted_genome_size > 6500000 ? "FLAG" : "PASS"' | \
    csvtk mutate2 -t -n coverage_qc \
        -e '$predicted_genome_size <= 0 ? "NA" : ($sum_len < 35 * $predicted_genome_size ? "FAIL" : ($sum_len < 40 * $predicted_genome_size ? "FLAG" : "PASS"))' | \
    csvtk mutate2 -t -n assembly_qc \
        -e '$predicted_genome_size <= 0 ? "NA" : ($assembly_length > 1.15 * $predicted_genome_size || $assembly_length < 0.85 * $predicted_genome_size ? "FAIL" : "PASS")' | \
    csvtk mutate2 -t -n contig_qc \
        -e '$contigs == 0 || $contigs > 30 ? "FAIL" : ($contigs > 10 ? "FLAG" : "PASS")' | \
    csvtk mutate -t -f species1 -n tmp_name1 -p '^(.*) \([0-9.]+%\)$' --na | \
    csvtk mutate -t -f species1 -n tmp_pct1 -p '\(([0-9.]+)%\)$' | \
    csvtk mutate -t -f species2 -n tmp_pct2 -p '\(([0-9.]+)%\)$' | \
    csvtk mutate -t -f species3 -n tmp_pct3 -p '\(([0-9.]+)%\)$' | \
    csvtk mutate -t -f species1 -n tmp_genus1 -p '^(\S+) ' | \
    csvtk mutate -t -f species2 -n tmp_genus2 -p '^(\S+) ' | \
    csvtk mutate -t -f species3 -n tmp_genus3 -p '^(\S+) ' | \
    csvtk mutate2 -t -n tmp_genus_pct -w 2 \
        -e '$tmp_pct1 + ($tmp_genus2 == $tmp_genus1 ? $tmp_pct2 : 0) + ($tmp_genus3 == $tmp_genus1 ? $tmp_pct3 : 0)' | \
    csvtk mutate2 -t -n species_qc \
        -e '$tmp_name1 == "NA" || $tmp_name1 == "" ? "FAIL" : ($tmp_genus_pct < 80 ? "FAIL" : ($tmp_pct1 < 80 ? "FLAG" : "PASS"))' \
    > ${OUTPUTDIR}/.temp_summary.tsv   # tmp_* helper columns are dropped when formatting below

echo 'Formatting summary report'

# select/reorder the raw summary columns (dropping lrge_qc) and give them informative headers.
# summary.tsv uses the same headers as summary.xlsx, minus the parenthesised details and with
# spaces replaced by underscores
SUMMARYCOLS_ORDER='file,reads,sum_len,min_len,avg_len,max_len,N50,species1,species2,species3,gsize_source,predicted_genome_size,mean_coverage,assembler,contigs,assembly_length,assembly_N50,coverage_qc,species_qc,assembly_qc,contig_qc'
SUMMARYCOLS_TSV='sequencing_ID,number_of_reads,sum_length,minimum_read_length,average_read_length,maximum_read_length,read_N50,kraken_top_species_call,kraken_species_call_2,kraken_species_call_3,predicted_genome_size_tool_used,predicted_genome_size,mean_depth,final_assembly_tool_used,number_of_contigs,assembly_length,assembly_N50,depth_QC,species_QC,assembly_QC,contig_QC'
SUMMARYCOLS_XLSX='sequencing ID,number of reads,sum length,minimum read length,average read length,maximum read length,read N50,kraken top species call,kraken species call 2,kraken species call 3,predicted genome size tool used,predicted genome size,mean depth (sum length / predicted genome size),final assembly tool used,number of contigs,assembly length,assembly N50,depth QC (mean coverage >40x = PASS; 35-40x = FLAG; <35x = FAIL),species QC (kraken top species call >80% = PASS; top 3 calls same genus <80% = FAIL; else FLAG),assembly QC (assembly length = predicted genome size +/- 15% PASS/FAIL),contig QC (>30 contigs = FAIL; 11-30 contigs = FLAG; <=10 contigs = PASS)'

csvtk cut -t -f "${SUMMARYCOLS_ORDER}" ${OUTPUTDIR}/.temp_summary.tsv | \
    csvtk rename -t -f "${SUMMARYCOLS_ORDER}" -n "${SUMMARYCOLS_TSV}" \
    > ${OUTPUTDIR}/summary.tsv

csvtk rename -t -f "${SUMMARYCOLS_TSV}" -n "${SUMMARYCOLS_XLSX}" ${OUTPUTDIR}/summary.tsv | \
    csvtk csv2xlsx -t -f -o ${OUTPUTDIR}/summary.xlsx

rm -f ${OUTPUTDIR}/.temp_summary.tsv

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
