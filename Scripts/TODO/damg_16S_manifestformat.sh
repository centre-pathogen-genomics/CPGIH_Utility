#!/bin/sh

i_flag=''
m_flag=''
o_flag=${PWD}
l_flag='Unspecified'
s_flag='false'
q_flag='false'
r_flag='Unspecified'
d_flag='false'
p_flag='skip'
a_flag='skip'
t_flag='8'
inputerror='false'

echo ''

while getopts ':hi:m:o:l:sqr:d:p:a:t:' opt
do
    case ${opt} in
	h)
	    echo -e 'Options:'
	    echo -e '\t -i Manifest File (TSV file with sample name and absolute paths to demultiplexed forward and reverse reads) [Mandatory]'
	    echo -e '\t -m Metadata File [Mandatory]'
	    echo -e '\t -o Output Directory [Optional]'
	    echo -e '\t\t If Unspecified: Will Output Files To Current Working Directory'
	    echo -e '\t -l Trim Length [Optional]'
	    echo -e '\t\t If Unspecified: Pipeline Will Pause And Wait for Prompt'
	    echo -e '\t\t This Is To Allow User To Evaluate Read Quality Metrics In demux.qzv'
	    echo -e '\t -s Deblur-Only Mode [Optional] [Default: False]'
            echo -e '\t\t Specifying This Option Will Skip Taxonomic Classification And Phylogenetic Placement Of sOTUs, Calculation Of Alpha And Beta Diversity, And ANCOM Analysis'
            echo -e '\t\t Useful When Samples Are Split Across Multiple Sequencing Runs'
	    echo -e '\t -q Calculation-Only Mode  [Optional] [Default: False]'
	    echo -e '\t\t Specifying This Option Will Skip Rarefaction, Calculation Of Alpha And Beta Diversity, And ANCOM Analysis'
	    echo -e '\t\t Taxonomic Classification And Phylogenetic Placement Of sOTUs Will Be Performed'
	    echo -e '\t\t Useful When Samples Are Split Across Multiple Sequencing Runs Or Downstream Analysis Will Be Performed By The User'
	    echo -e '\t -r Rarefaction Depth [Optional]'
	    echo -e '\t -d Run Decontam [Optional] [Default: False]'
	    echo -e '\t\t If Desired: Specify Name Of Column In Metadata Sheet Describing Whether Each Sample Is A Negative Control'
	    echo -e '\t\t Column Must Be Boolean With Values TRUE Or FALSE'
	    echo -e '\t -p Metadata Variable(s) For PERMANOVA [Optional] [Default: Skip]'
	    echo -e '\t\t For Multiple Variables Use A Space Separated List Within Double Quotes'
	    echo -e '\t\t eg. "meta1 meta2 meta3"'
	    echo -e '\t -a Metadata Variable(s) For ANCOM [Optional] [Default: Skip]'
	    echo -e '\t\t eg. "meta1 meta2 meta3"'
	    echo -e '\t\t For Multiple Variables Use A Space Separated List Within Double Quotes'
	    echo -e '\t -t Threads Or Parallel Jobs To Use When Possible [Optional] [Default: 8]'
	    echo ''
	    exit 0
	    ;;
	i) i_flag=$OPTARG ;;
	m) m_flag=$OPTARG ;;
	o) o_flag=$OPTARG ;;
	l) l_flag=$OPTARG ;;
	r) r_flag=$OPTARG ;;
	d) d_flag=$OPTARG ;;
	p) p_flag=$OPTARG ;;
	a) a_flag=$OPTARG ;;
	s) s_flag='true' ;;
	q) q_flag='true' ;;
	t) t_flag=$OPTARG ;;
	\?) echo -e '\tUsage: damg_16S_manifestformat.sh -i PathToManifestFile -m MetadataFile\n\tOR\n\tHelp and Optional Arguments: damg_16S_manifestformat.sh -h\n' >&2
	    exit 1
	    ;;
	:) echo -e '\tError: Use -h for full options list\n'
	   exit 1
    esac
done

# ENSURE MANIFEST AND METADATA FILES ARE SPECIFIED
# IF FALSE: PRINT ERROR AND EXIT
# IF TRUE: PRINT THEIR NAMES, THREADS SPECIFIED, AND OUTPUT DIRECTORY PATH
if [ "$i_flag" = '' ] || [ "$m_flag" = '' ]
then
    echo 'ERROR'
    echo 'Usage: damg_16S_manifestformat.sh -i ManifestFile -m MetadataFile'
    echo 'Use -h for full options list'
    echo ''
    exit 1
else
    echo 'Manifest file: '$i_flag
    echo 'Metadata File: '$m_flag
    echo 'Threads or Parallel Jobs To Use: '$t_flag
    if [ "$o_flag" = ${PWD} ]
    then
		echo 'Writing Output To Current Working Directory'
    else
		echo 'Output Directory: '$o_flag
    fi
fi

# PRINT TRIM LENGTH OR ASK FOR USER INPUT
if [ $l_flag = 'Unspecified' ]
then
    echo 'Trim Length: Unspecified - User Input Will Be Required'
else
    echo 'Trim Length: '$l_flag'bp'
fi

# PRINT WHETHER DECONTAM STEP IS REQUESTED AND NAME OF COLUMN IN METADATA
if [ $d_flag = 'false' ]
then
    echo 'Skipping Decontam'
else
    echo 'Decontam metadata column: '$d_flag
fi

# IF DEBLUR-ONLY RUN MODE WAS NOT SPECIFIED
# PRINT RAREFACTION DEPTH OR ASK FOR USER INPUT
# PRINT WHETHER PERMOVA AND ANCOM WILL RUN, AND ON WHICH VARIABLES
if [ $s_flag = 'false' ]
then
    if [ $r_flag = 'Unspecified' ]
    then
		echo -e 'Rarefaction Depth: Unspecified - User Input May Be Required'
    else
		echo -e 'Rarefaction Depth: '$r_flag
    fi

    if [[ $p_flag = 'skip' ]]
    then
		echo -e 'Metadata Variable(s) for PERMANOVA: Unspecified - Skipping'
    else
		echo -e 'Metadata Variable(s) for PERMANOVA: '$p_flag
    fi

    if [[ $a_flag = 'skip' ]]
    then
		echo -e 'Metadata Variable(s) for ANCOM: Unspecified - Skipping'
    else
		echo -e 'Metadata Variable(s) for ANCOM: '$a_flag
    fi
elif [ $s_flag = 'true' ]
then
    echo 'Deblur-Only Run Mode Specified (Pipeline Will Finish After sOTU generation)'
else
    echo 'Deblur-Only Run Mode Incorrectly Specified'
    inputerror='true'
fi

# PRINT WHETHER CALCULATION-ONLY RUN MODE WAS SPECIFIED
if [ $q_flag = 'false' ]
then
    :
elif [ $q_flag = 'true' ]
then
    echo 'Calculation-Only Run Mode Specified (Pipeline Will Skip Rarefaction, Diversity Calculations, Statistics, and ANCOM)'
else
    echo 'Calculation-Only Run Mode Incorrectly Specified'
    inputerror='true'
fi  

# CONFIRM THAT MANIFEST AND METADATA FILES EXIST
if [ ! -f $i_flag ] || [ ! -f $m_flag ]
then
    if [ ! -f $i_flag ]
    then
		echo 'Manifest File Not Found'
		inputerror='true'
    fi

    if [ ! -f $m_flag ]
    then
		echo 'Metadata File Not Found'
		inputerror='true'
    fi
fi

# IF OUTPUT DIRECTORY DOESNT EXIT, THEN CREATE
# IF IT IS CURRENT WORKING DIRECTORY, DO NOTHING
# IF OUTPUT DIRECTORY ALREADY EXISTS, WARN THAT CONTENTS WILL BE OVERWRITTEN
if [ ! -d $o_flag ]
then
    mkdir -p $o_flag
elif [ "$o_flag" = ${PWD} ]
then
    :
else
    echo 'Output Directory '$o_flag' Already Exists: Contents Will Be Overwritten'
fi

# CHECK WHETHER MANIFEST HEADERS ARE CORRECT - QIIME2 IS VERY FUSSY ABOUT THESE
# IF THEY ARE INCORRECT, FIX
if $(head -n 1 $i_flag | grep -qP '^sample-id\tforward-absolute-filepath\treverse-absolute-filepath$')
then
	:
else
	echo 'Manifest input headers are incorrect. Fixing.'
	sed -i '1s/^/sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n/' $i_flag
fi

# CHECK THE VALUE IN [1,1] OF METADATA FILE
# QIIME2 IS FUSSY ABOUT WHAT VALUES IT LIKES
# CHANGE TO 'sample_name' TO MAKE SURE IT IS HAPPY
awk 'BEGIN{FS=OFS="\t"} NR==1{$1="sample_name"}1' $m_flag > tmpfile && mv -f tmpfile $m_flag

# CONFIRM EXISTENCE OF ALL FASTQ FILES SPECIFIED IN MANIFEST
sed '1d' $i_flag > manifest.temp

while read sample read1 read2
do
	if [ ! -f $read1 ]
	then
		echo 'File '$read1' does not exist'
		inputerror='true'
	fi

	if [ ! -f $read2 ]
	then
		echo 'File '$read2' does not exist'
		inputerror='true'
	fi
done < manifest.temp

rm -f manifest.temp
  
 # PRINT EMPTY LINE FOR FORMATTING 
echo ''

# IF ANY ABOVE TESTS THREW ERRORS, STOP PIPELINE HERE
if [ $inputerror = 'true' ]
then
    exit 1
fi

##### CHECKS COMPLETE, STARTING PIPELINE

# IMPORT DEMULTIPLEXED FASTQ READS AS QIIME2 ARTEFACT
qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
  	--input-path $i_flag \
  	--output-path demux.qza \
  	--input-format PairedEndFastqManifestPhred33V2

# SUMMARISE READS PER SAMPLE ETC.
qiime demux summarize \
	--i-data demux.qza \
	--o-visualization demux.qzv

# QUALITY CONTROL OF FASTQ READS
qiime quality-filter q-score \
	--i-demux demux.qza \
	--o-filtered-sequences demux-filtered.qza \
	--o-filter-stats demux-filterstats.qza

# SUMMARISE QC STATS
qiime metadata tabulate \
	--m-input-file demux-filterstats.qza \
	--o-visualization demux-filterstats.qzv

# RUN DEBLUR TO DEFINE SUB-OTUS
# ALL READS NEED TO BE TRIMMED TO THE SAME LENGTH
# WILL TAKE VALUE FROM COMMAND LINE IF SPECIFIED
# OR ELSE PAUSE AND ASK FOR KEYBOARD INPUT 
if [ $l_flag = 'Unspecified' ]
then
    read -p 'Trim Length: ' trimlength_var
    qiime deblur denoise-16S \
		--i-demultiplexed-seqs demux-filtered.qza \
		--p-trim-length $trimlength_var \
		--p-sample-stats \
		--p-jobs-to-start $t_flag \
		--p-no-hashed-feature-ids \
		--o-table feature-table.qza \
		--o-representative-sequences rep-seqs.qza \
		--o-stats deblur-stats.qza
else
    qiime deblur denoise-16S \
		--i-demultiplexed-seqs demux-filtered.qza \
		--p-trim-length $l_flag \
		--p-sample-stats \
		--p-jobs-to-start $t_flag \
		--p-no-hashed-feature-ids \
		--o-table feature-table.qza \
		--o-representative-sequences rep-seqs.qza \
		--o-stats deblur-stats.qza
fi

# RUN DEBLUR TO REMOVE CONTAMINANT OTUS IF SPECIFIED
if [ $d_flag != 'false' ]
then
	mv feature-table.qza feature-table-predecontam.qza

	Rscript "$(dirname "$(realpath "$0")""/decontam.R feature-table-predecontam.qza feature-table-decontam.biom $m_flag $d_flag

	qiime tools import \
		--input-path feature-table-decontam.biom \
		--type 'FeatureTable[Frequency]' \
		--input-format BIOMV100Format \
		--output-path feature-table.qza

	rm feature-table-decontam.biom
fi

# END PIPELINE HERE IF IF DEBLUR-ONLY RUN MODE WAS SPECIFIED
if [ $s_flag = 'true' ]
then
	# MOVE OUTPUT FILES TO OUTPUT DIRECTORY IF SPECIFIED
	if [ ! $o_flag = $PWD ]
	then
		mv deblur.log deblur-stats.qza $o_flag
		mv demux-filtered.qza demux-filterstats.qza demux-filterstats.qzv demux.qza demux.qzv $o_flag
		mv feature-table.qza rep-seqs.qza $o_flag
		rm -f raw-seqs.qza

		if [ $d_flag != 'false' ]
		then
			mv contaminant_otus.csv feature-table-predecontam.qza $o_fla
		fi
	fi

    echo 'Skipping Analysis Steps As Instructed'
    echo 'Pipeline Complete!'
    echo ''
    exit 0
fi

# SUMMARISE FEATURES PER SAMPLE ETC
qiime feature-table summarize \
	--i-table feature-table.qza \
	--m-sample-metadata-file $m_flag \
	--o-visualization feature-table.qzv

# MAKE A VISUALISATION OF OTUS SEQUENCES 
qiime feature-table tabulate-seqs \
	--i-data rep-seqs.qza \
	--o-visualization rep-seqs.qzv

# CONSTRUCT OTU PHYLOGENY USING SEPP METHOD
qiime fragment-insertion sepp \
	--i-representative-sequences rep-seqs.qza \
	--i-reference-database /home/cwwalsh/Scripts/DAMG/sepp-refs-gg-13-8.qza \
	--o-tree insertion-tree.qza \
	--o-placements insertion-placements.qza

# ASSIGN TAXONOMY TO OTUS USING PRETRAINED V4 CLASSIFIER MODEL
qiime feature-classifier classify-sklearn \
	--i-reads rep-seqs.qza \
	--i-classifier /home/cwwalsh/Scripts/DAMG/gg-13-8-99-515-806-nb-classifier.qza \
	--o-classification taxonomy.qza

# MAKE VISUALISTION OF OTU TAXONOMY PER SEQUENCE
qiime metadata tabulate \
	--m-input-file taxonomy.qza \
	--o-visualization taxonomy.qzv

# EXPORT TAXONOMY TABLE AS TSV
qiime tools export \
	--input-path taxonomy.qza \
	--output-path .

# EXPORT TREE AS NEWICK
qiime tools export \
	--input-path insertion-tree.qza \
	--output-path .

# MODIFY TREE DELIMITERS TO BE APE/PHYLOSEQ COMPLIANT
sed -i 's/; /| /' tree.nwk

# EXPORT OTU TABLE AS BIOM
qiime tools export \
	--input-path feature-table.qza \
	--output-path .

# COMVERT OTU TABLE BIOM TO JSON 
biom convert \
	-i feature-table.biom \
	-o feature-table_json.biom \
	--table-type="OTU table" \
	--to-json

# MAKE VISUALISATION OF TAXONOMIC BARPLOTS FOR EACH SAMPLE
qiime taxa barplot \
	--i-table feature-table.qza \
	--i-taxonomy taxonomy.qza \
	--m-metadata-file $m_flag \
	--o-visualization taxa-bar-plots.qzv

# END PIPELINE HERE IF CALCULATION-ONLY RUN MODE WAS SPECIFIED
if [ $q_flag = 'true' ]
then
	# MOVE OUTPUT FILES TO OUTPUT DIRECTORY IF SPECIFIED
    if [ ! $o_flag = $PWD ]
    then
		mv deblur.log deblur-stats.qza $o_flag
        mv demux-filtered.qza demux-filterstats.qza demux-filterstats.qzv demux.qza demux.qzv $o_flag
        mv feature-table.qza feature-table.qzv rep-seqs.qza $o_flag
        mv rep-seqs.qzv insertion-tree.qza insertion-placements.qza taxonomy.qza taxonomy.qzv $o_flag
        mv alpha-rarefaction-*.qzv taxa-bar-plots.qzv $o_flag

		# MOVE DECONTAM FILES TO OUTPUT DIRECTORY IF SPECIFIED
		if [ $d_flag != 'false' ]
		then
	    	mv contaminant_otus.csv feature-table-predecontam.qza $o_fla
		fi
	fi
		
    echo 'Skipping Analysis Steps As Instructed'
    echo 'Pipeline Complete!'
    echo ''
    exit 0
fi

# PERFORM RAREFACTION AND DIVERSITY CACULATIONS
errcode=1
while [ $errcode -ne 0 ]
do
	# IF DEPTH IS UNSPECIFIED, ASK FOR KEYBOARD INPUT, CREATE RAREFACTION CURVE
	# CAN TRY MULTIPLE DEPTHS WITH KEYBOARD INPUT TO FIND A SUITABLE ONE
	# BEFORE MOVING TO DIVERSITY CALCULATIONS
    if [ $r_flag = 'Unspecified' ]
    then
		
		read -p 'Rarefaction depth(s) for testing (can provide multiple depths separated by a space): ' rarefactiondepths_var
	
	for depth in $rarefactiondepths_var
	do
	    
		qiime diversity alpha-rarefaction \
			--i-table feature-table.qza \
			--i-phylogeny insertion-tree.qza \
			--p-max-depth $depth \
			--m-metadata-file $m_flag \
			--o-visualization alpha-rarefaction-"$depth".qzv
	
	done
	
	errcode=$?
	
	# IF RAREFACTION DEPTH INPUT IF NOT ALLOWABLE, ASK USER TO RESPECIFY
	if [ $errcode -ne 0 ]
	then
	    echo 'Inappropriate Rarefaction Depth Specified'
	    echo 'Please Retry'
	    errcode=1
	# ASK FOR DEPTH FOR DIVERSITY CALCULATIONS
	else
	    echo 'Rarefaction testing complete!'
	    read -p 'Rarefaction depth for diversity analysis (entering no value will skip this step): ' diversityrarefactiondepth_var

		# PERFORM DIVERSITY CALCULATIONS
		# IF DEPTH IS SPECIFIED
		# SKIP IF NOT
	    if [ ! $diversityrarefactiondepth_var = '' ]
	    then
		qiime diversity core-metrics-phylogenetic \
			--i-table feature-table.qza \
			--i-phylogeny insertion-tree.qza \
			--p-sampling-depth $diversityrarefactiondepth_var \
			--m-metadata-file $m_flag \
			--p-n-jobs-or-threads $t_flag \
			--output-dir CoreMetricsPhylogenetic
	    fi
	fi
    else

	# MAKE RAREFACTION CURVE TO SPECIFIED DEPTH
	qiime diversity alpha-rarefaction \
		--i-table feature-table.qza \
		--i-phylogeny insertion-tree.qza \
		--p-max-depth $r_flag \
		--m-metadata-file $m_flag \
		--o-visualization alpha-rarefaction-"$r_flag".qzv

	# CACULATE DIVERSITY METRICS AT SPECIFIED DEPTH
	qiime diversity core-metrics-phylogenetic \
		--i-table feature-table.qza \
		--i-phylogeny insertion-tree.qza \
		--p-sampling-depth $r_flag \
		--m-metadata-file $m_flag \
		--p-n-jobs-or-threads $t_flag \
		--output-dir CoreMetricsPhylogenetic

	errcode=$?
    fi
done

# ALPHA DIVERSITY STATISTICS 
for i in $(ls CoreMetricsPhylogenetic/*_vector.qza | sed 's/_vector\.qza//')
do
    qiime diversity alpha-group-significance \
		--i-alpha-diversity "$i"_vector.qza \
		--m-metadata-file $m_flag \
		--o-visualization "$i"_groupsig.qzv
done

# BETA DIVERSITY STATISTICS IF SPECIFIED
if [[ $p_flag = 'skip' ]]
then
	echo 'Skipping PERMANOVA As Instructed'
else
    for metric in $(ls CoreMetricsPhylogenetic/*_distance_matrix.qza | sed 's/_distance_matrix.qza//')
    do
		for variable in $p_flag
		do
	    	qiime diversity beta-group-significance \
				--i-distance-matrix "$metric"_distance_matrix.qza \
				--m-metadata-file $m_flag \
				--m-metadata-column $variable \
				--p-pairwise \
				--o-visualization "$metric"_groupsig_"$variable".qzv
		done
    done
fi

# ANCOM ANALYSIS IF SPECIFIED
if [[ $a_flag = 'skip' ]]
then
    echo 'Skipping ANCOM As Instructed'
else
	qiime composition add-pseudocount \
		--i-table feature-table.qza \
		--o-composition-table comp-feature-table.qza

    for variable in $a_flag
    do
		qiime composition ancom \
			--i-table comp-feature-table.qza \
			--m-metadata-file $m_flag \
			--m-metadata-column $variable \
			--o-visualization ancom-"$variable".qzv
    done
fi

# EXPORT DISTANCE MATRICES 
for metric in $(ls CoreMetricsPhylogenetic/*_distance_matrix.qza | sed 's/_distance_matrix.qza// ; s/CoreMetricsPhylogenetic\///')
do
    qiime tools export \
        --input-path CoreMetricsPhylogenetic/"$metric"_distance_matrix.qza \
        --output-path CoreMetricsPhylogenetic/ \
    
    mv CoreMetricsPhylogenetic/distance-matrix.tsv CoreMetricsPhylogenetic/"$metric"_distance_matrix.tsv
done

# MOVE OUTPUT FILES TO OUTPUT DIRECTORY IF SPECIFIED
if [ ! $o_flag = $PWD ]
then
	mv deblur.log deblur-stats.qza $o_flag
	mv demux-filtered.qza demux-filterstats.qza demux-filterstats.qzv demux.qza demux.qzv $o_flag
	mv feature-table.qza feature-table.qzv rep-seqs.qza $o_flag
	mv rep-seqs.qzv insertion-tree.qza insertion-placements.qza taxonomy.qza taxonomy.qzv $o_flag
	mv alpha-rarefaction-*.qzv CoreMetricsPhylogenetic/ taxa-bar-plots.qzv comp-feature-table.qza $o_flag

	if [ $d_flag != 'false' ]
		then
			mv contaminant_otus.csv feature-table-predecontam.qza $o_flag
		fi
	
	if [ ! $a_flag = '' ] && [ ! $ancomvariables_var = '' ]
	then
		mv ancom-*.qzv $o_flag
	fi

	mv taxonomy.tsv tree.nwk feature-table.biom feature-table_json.biom $o_flag
	
	rm -f raw-seqs.qza 
fi
	
echo 'Pipeline Complete!'
echo ''
