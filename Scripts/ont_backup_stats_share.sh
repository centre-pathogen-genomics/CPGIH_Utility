#!/bin/sh

# USAGE: /home/cwwalsh/Scripts/DAMG/ont_backup_stats_share.sh runID renamingTSV /optional/path/to/share_template.txt

set -e
RUNID=$1
RENAMINGTSV=$2
SHARETEMPLATE=$3

INPUTDIR=/home/damg/for_transfer/${RUNID}
OUTPUTDIR=/home/damg/data/${RUNID}

# CONFIRM THAT INPUT DIR EXISTS
if [ ! -d ${INPUTDIR} ]
then

    echo "Input Directory Does Not Exist"
    exit 1

else

    echo "Input Directory:" ${INPUTDIR}

fi

# CONFIRM THAT OUTPUT DIR DOES NOT EXIST
if [ -d ${OUTPUTDIR} ]
then

    echo "Output Directory Already Exists"
    exit 1

fi

# CONFIRM THAT SHARING TEMPLATE FILE EXISTS IF PROVIDED
if [ -z ${3} ]
then
    
    echo "Sharing File Template Not Specified. Will Not Generate Share Link"

else
    if [ ! -f ${SHARETEMPLATE} ]
    then

        echo "Sharing File Template Does Not Exist"
        exit 1

    fi

fi

# CONFIRM INPUT FASTQ DIR EXISTS
if [ ! -d ${INPUTDIR}/fastq_pass ] && [ ! -d ${INPUTDIR}/fastq ] && [ ! -d ${INPUTDIR}/fastq-pass ]
then

    echo "Input FASTQ Directory Does Not Exist"
    exit 1

fi

# CONFIRM INPUT POD5 DIR EXISTS
if [ ! -d ${INPUTDIR}/pod5_pass ] && [ ! -d ${INPUTDIR}/pod5 ]
then

    echo "Input POD5 Directory Does Not Exist"
    exit 1

fi

# CONFIRM INPUT REPORT DIR EXISTS
if [ ! -d ${INPUTDIR}/reports ] && [ ! -d ${INPUTDIR}/other_reports ]
then

    echo "Input Reports Directory Does Not Exist"
    exit 1

fi

# IF WRITE PERMISSIONS ARE SET CORRECTNLY ON INPUT DIR
# MAKE OUTPUT DIR
# MOVE FASTQ, POD5, REPORT, AND ANALYSIS DIRS TO OUTPUT DIR
if [ -w ${INPUTDIR} ]
then

    echo "Ouput Directory:" ${OUTPUTDIR}
    mkdir ${OUTPUTDIR}
    
    if [ -d ${INPUTDIR}/fastq_pass ]
    then

        mv ${INPUTDIR}/fastq_pass ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/fastq ]
    then

        mv ${INPUTDIR}/fastq ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/fastq-pass ]
    then

        mv ${INPUTDIR}/fastq-pass ${OUTPUTDIR}

    fi
    
    if [ -d ${INPUTDIR}/fast5_pass ]
    then

        mv ${INPUTDIR}/fast5_pass ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/pod5_pass ]
    then

        mv ${INPUTDIR}/pod5_pass ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/fast5 ]
    then

        mv ${INPUTDIR}/fast5 ${OUTPUTDIR}

    fi
    
    if [ -d ${INPUTDIR}/pod5 ]
    then

        mv ${INPUTDIR}/pod5 ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/reports ]
    then

        mv ${INPUTDIR}/reports ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/other_reports ]
    then

        mv ${INPUTDIR}/other_reports ${OUTPUTDIR}

    fi

    if [ -d ${INPUTDIR}/analysis ]
    then

        mv ${INPUTDIR}/analysis ${OUTPUTDIR}

    fi
else

    echo "Permission Denied"
    exit 1

fi

# COPY RENAMING INFO FILE TO OUTPUT DIR
cp ${RENAMINGTSV} ${OUTPUTDIR}/renaming.tsv

# SET WORKING DIR TO OUTPUT FASTQ OR FASTQ_PASS DIR 
if [ -d ${OUTPUTDIR}/fastq_pass/ ]
then

    cd ${OUTPUTDIR}/fastq_pass/

elif [ -d ${OUTPUTDIR}/fastq/ ]
then

    cd ${OUTPUTDIR}/fastq/

fi

# CONCATENATE FASTQ IF NEEDED
# RENAME
# DELETE ORIGINAL BARCODE DIRECTORY
SAMPLECOUNT=$(ls -d * | wc -l)
echo ${SAMPLECOUNT} "Samples Found"
    
while IFS= read -r i j 
do

    FILECOUNT=$(ls -1 ${i}/ | wc -l)

    if [ ${FILECOUNT} == 1 ]
    then

        echo "One FASTQ For Sample" ${i} "No Concatentation Needed. Just Renaming."
        mv ${i}/*${i}*.fastq.gz ${j}.fastq.gz
	rmdir ${i}

    else

        echo "Concatenating" ${FILECOUNT} "FASTQ From" ${i} "And Renaming."
        pigz -cd -p 10 ${i}/*.fastq.gz | pigz -p 10 > ${j}.fastq.gz
	rm -rf ${i}

    fi

done < <(cat ../renaming.tsv; echo)

# CREATE SEQ STATS
echo "Generating Read Length and Quality Statistics"
seqkit stats -abT *.fastq.gz | \
    cut -f 1,4,5,13,14,15 | \
        csvtk pretty -t > ../seqkit_stats.tsv

# SHARE ALL FILES FROM RUN IF SHARING TEMPLATE FILE PROVIDED
if [ ! -z ${3} ]
then

    ls *.fastq.gz | sed 's/.fastq.gz$//' > names.tsv
    mdu share --input_file names.tsv --source . -t ${SHARETEMPLATE}
    rm -f names.tsv

fi

# PRINT SEQ STATS TO SCREEN
cat ../seqkit_stats.tsv

# REMOVE EMPTY INPUT DIRECTORY
rmdir ${INPUTDIR}
