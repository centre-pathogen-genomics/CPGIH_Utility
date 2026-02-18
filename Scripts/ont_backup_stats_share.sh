#!/bin/sh

# USAGE: ont_backup_stats_share.sh runID renamingTSV /optional/path/to/share_template.txt

set -e
RUNID=$1
RENAMINGTSV=$2
SHARETEMPLATE=$3

INPUTDIRECTORY=/home/damg/for_transfer/"$RUNID"
OUTPUTDIRECTORY=/home/damg/data/"$RUNID"

# CONFIRM THAT INPUT DIR EXISTS
if [ ! -d "$INPUTDIRECTORY" ]
then

    echo "Input Directory Does Not Exist"
    exit 1

else

    echo "Input Directory:" "$INPUTDIRECTORY"

fi

# CONFIRM THAT OUTPUT DIRECTORY DOES NOT EXIST
if [ -d "$OUTPUTDIRECTORY" ]
then

    echo "Output Directory Already Exists"
    exit 1

fi

# CONFIRM THAT SHARING TEMPLATE FILE EXISTS IF PROVIDED
if [ -z "$3" ]
then
    
    echo "Sharing File Template Not Specified. Will Not Generate Share Link"

else
    if [ ! -f "$SHARETEMPLATE" ]
    then

        echo "Sharing File Template Does Not Exist"
        exit 1

    fi

fi

# CONFIRM INPUT FASTQ DIRECTORY EXISTS
if [ ! -d "$INPUTDIRECTORY"/fastq_pass ] && [ ! -d "$INPUTDIRECTORY"/fastq ] && [ ! -d "$INPUTDIRECTORY"/fastq-pass ]
then

    echo "Input FASTQ Directory Does Not Exist"
    exit 1

fi

# CONFIRM INPUT POD5 DIRECTORY EXISTS
if [ ! -d "$INPUTDIRECTORY"/pod5_pass ] && [ ! -d "$INPUTDIRECTORY"/pod5 ]
then

    echo "Input POD5 Directory Does Not Exist"
    exit 1

fi

# CONFIRM INPUT REPORT DIRECTORY EXISTS
if [ ! -d "$INPUTDIRECTORY"/reports ] && [ ! -d "$INPUTDIRECTORY"/other_reports ]
then

    echo "Input Reports Directory Does Not Exist"
    exit 1

fi

# IF WRITE PERMISSIONS ARE SET CORRECTNLY ON INPUT DIRECTORY
# MAKE OUTPUT DIRECTORY AND COPY RENAMING FILE 
# OTHERWISE, ERROR AND EXIT
if [ -w "$INPUTDIRECTORY" ]
then

    echo "Ouput Directory:" "$OUTPUTDIRECTORY"
    mkdir "$OUTPUTDIRECTORY"

    cp "$RENAMINGTSV" "$OUTPUTDIRECTORY"/renaming.tsv

else

    echo "Permission Denied"
    exit 1

fi

# MOVE POD5, REPORT, AND ANALYSIS DIRECTORIES TO OUTPUT DIRECTORY
if [ -d "$INPUTDIRECTORY"/fast5_pass ]
then

    mv "$INPUTDIRECTORY"/fast5_pass "$OUTPUTDIRECTORY"

fi

if [ -d "$INPUTDIRECTORY"/pod5_pass ]
then

    mv "$INPUTDIRECTORY"/pod5_pass "$OUTPUTDIRECTORY"

fi

if [ -d "$INPUTDIRECTORY"/fast5 ]
then

    mv "$INPUTDIRECTORY"/fast5 "$OUTPUTDIRECTORY"

fi

if [ -d "$INPUTDIRECTORY"/pod5 ]
then

    mv "$INPUTDIRECTORY"/pod5 "$OUTPUTDIRECTORY"

fi

if [ -d "$INPUTDIRECTORY"/reports ]
then

    mv "$INPUTDIRECTORY"/reports "$OUTPUTDIRECTORY"

fi

if [ -d "$INPUTDIRECTORY"/other_reports ]
then

    mv "$INPUTDIRECTORY"/other_reports "$OUTPUTDIRECTORY"

fi

if [ -d "$INPUTDIRECTORY"/analysis ]
then

    mv "$INPUTDIRECTORY"/analysis "$OUTPUTDIRECTORY"

fi

# ONLY MOVE FASTQ DATA FOR SAMPLES NAMED IN RENAMINGTSV
# AND CHANGE WORKING DIRECTORY TO FASTQ DIRECTORY
if [ -d "$INPUTDIRECTORY"/fastq_pass ]
then

    mkdir "$OUTPUTDIRECTORY"/fastq_pass

    while IFS=$'\t' read -r i j || [[ -n ${i} ]] 
    do

        mv "$INPUTDIRECTORY"/fastq_pass/"$i" "$OUTPUTDIRECTORY"/fastq_pass/
    
    done < "$RENAMINGTSV"

    cd "$OUTPUTDIRECTORY"/fastq_pass/

fi

if [ -d "$INPUTDIRECTORY"/fastq ]
then

    mkdir "$OUTPUTDIRECTORY"/fastq

    while IFS=$'\t' read -r i j || [[ -n ${i} ]] 
    do

        mv "$INPUTDIRECTORY"/fastq/"$i" "$OUTPUTDIRECTORY"/fastq/
    
    done < "$RENAMINGTSV"

    cd "$OUTPUTDIRECTORY"/fastq/

fi

if [ -d "$INPUTDIRECTORY"/fastq-pass ]
then

    mkdir "$OUTPUTDIRECTORY"/fastq-pass

    while IFS=$'\t' read -r i j || [[ -n ${i} ]] 
    do

        mv "$INPUTDIRECTORY"/fastq-pass/"$i" "$OUTPUTDIRECTORY"/fastq-pass/
    
    done < "$RENAMINGTSV"

    cd "$OUTPUTDIRECTORY"/fastq-pass/

fi

# CONCATENATE FASTQ IF NEEDED
# RENAME
# DELETE ORIGINAL BARCODE DIRECTORY
SAMPLECOUNT=$(ls -d * | wc -l)
echo "$SAMPLECOUNT" "Samples Found"
    
while IFS=$'\t' read -r i j || [[ -n ${i} ]] 
do

    FILECOUNT=$(ls -1 ${i}/ | wc -l)

    if [ $FILECOUNT = 1 ]
    then

        echo "One FASTQ For Sample" ${i} "No Concatentation Needed. Just Renaming."
        mv ${i}/*${i}*.fastq.gz ${j}.fastq.gz
	    rmdir ${i}

    else

        echo "Concatenating" "$FILECOUNT" "FASTQ From" ${i} "And Renaming."
        pigz -cd -p 10 ${i}/*.fastq.gz | pigz -p 10 > ${j}.fastq.gz
	    rm -rf ${i}

    fi

done < ../renaming.tsv

# CREATE SEQ STATS
echo "Generating Read Length and Quality Statistics"
/home/cwwalsh/Software/seqkit stats -abT *.fastq.gz | \
    cut -f 1,4,5,13,14,15 | \
        /home/cwwalsh/Software/csvtk pretty -t > ../seqkit_stats.tsv

# SHARE ALL FASTQ FILES FROM RUN IF SHARING TEMPLATE FILE PROVIDED
if [ ! -z "$3" ]
then

    ls *.fastq.gz | sed 's/.fastq.gz$//' > names.tsv
    /home/mdu/bin/mdu share --input_file names.tsv --source . -t "$SHARETEMPLATE"
    rm -f names.tsv

fi

# TREATING BAM FILES THE SAME AS FASTQ FILES
# BACKING UP AND RENAMING IF LISTED IN RENAMINGTSV
if [ -d "$INPUTDIRECTORY"/bam_pass ]
then

    mkdir "$OUTPUTDIRECTORY"/bam_pass

    while IFS=$'\t' read -r i j || [[ -n ${i} ]] 
    do
      
        mv "$INPUTDIRECTORY"/bam_pass/"$i"/*.bam "$OUTPUTDIRECTORY"/bam_pass/"$j".bam
    
    done < "$RENAMINGTSV"

fi

# PRINT SEQ STATS TO SCREEN
cat ../seqkit_stats.tsv

# REMOVE EMPTY INPUT DIRECTORY
rmdir "$INPUTDIRECTORY"
