#!/bin/bash

# USAGE: bash ont_combine_fastq.sh inputdirectory outputdirectory renamingtsv outputformat

set -e
INDIR=$1
OUTDIR=$2
RENAMING=$3

# CONFIRM THAT INDIR EXISTS
if [ ! -d ${INDIR} ]
then

    echo "Input Directory Does Not Exist"
    exit 1

else

    echo "Input Directory:" ${INDIR}

fi

# CONFIRM THAT OUTPUT DIR DOES NOT EXIST
if [ -d ${OUTDIR} ]
then

    echo "Output Directory Already Exists"
    exit 1
else

    echo "Output Directory:" ${OUTDIR}

fi

# CONFIRM THAT RENAMING FILE EXISTS
if [ ! -f ${RENAMING} ]
then

    echo "Renaming File Does Not Exist"
    exit 1

else

    echo "Renaming File:" ${RENAMING}

fi

# CONFIRM THAT OUTPUT FORMAT IS SPECIFIED AND VALID
if [ -z $4 ]
then

    OUTFORMAT=$4

    if [[ ${OUTFORMAT} == "file" || ${OUTFORMAT} == "subdir" ]]
    then

        echo "Specified Output Format:" ${OUTFORMAT}

    else

        echo "Specified Output Format Is Invalid, Defaulting to \"file\""
        OUTFORMAT="file"

    fi

else

    echo "Output Format Unspecified. Defaulting to \"file\""
    OUTFORMAT="file"

fi

# MAKE OUTDIR
mkdir -p ${OUTDIR}

# CONCATENATE FASTQ AND RENAME 
while IFS='\t' read -r i j || [[ -n "$i" ]]
do

    FILECOUNT=$( ls -1 ${INDIR}/${i}/ | wc -l )

    if [ ${FILECOUNT} -eq 1 ]
    then

        echo "One FASTQ For Sample" ${i} "No Concatentation Needed. Just Renaming (" ${j} ")"

        if [ ${OUTFORMAT} == "file" ]
        then

            mv ${INDIR}/${i}/*.fastq.gz ${OUTDIR}/${j}.fastq.gz

        else

            mkdir -p ${OUTDIR}/${j}/
            mv ${INDIR}/${i}/*.fastq.gz ${OUTDIR}/${j}/${j}.fastq.gz

        fi

    else

        echo "Concatenating" ${FILECOUNT} "FASTQ From" ${i} "And Renaming (" ${j} ")"

        if [[ ${OUTFORMAT} == "file" ]]
        then

            pigz -cd -p 10 ${INDIR}/${i}/*.fastq.gz | pigz -p 10 > ${OUTDIR}/${j}.fastq.gz

        else

            mkdir -p ${OUTDIR}/${j}/
            pigz -cd -p 10 ${INDIR}/${i}/*.fastq.gz | pigz -p 10 > ${OUTDIR}/${j}/${j}.fastq.gz

        fi

    fi

done < ${RENAMING}

