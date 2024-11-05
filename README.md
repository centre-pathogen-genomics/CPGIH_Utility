# Centre for Pathogen Genomics: Innovation Hub Utility Scripts

This repository contains utility scripts for random stuff we do in CPG-IH.  

## INSTALLATION AND SETUP
First you need to clone this repository to your local machine/server. You can do this to any directory, but for the purposes of this documentation we will do it in a new directory called `Tools/` in your home directory (`~`).  

```bash
mkdir ~/Tools/
cd ~/Tools/
git clone https://github.com/centre-pathogen-genomics/CPGIH_Utility.git 
```

If you are using the MDU servers is it easiest to just load my conda env - it should set all the paths correctly.  

```bash
conda activate /home/cwwalsh/miniconda3/envs/cpgih_utility
```

Alternatively you can make your own conda environment. Following these steps should give you one that does everything - you will need to install your own databases for `kraken2` and `emu` or modify the scripts to use existing ones.  

```bash
conda create -n cpgih_utility -y
conda activate cpgih_utility
conda install -c bioconda kraken2 shovill seqkit csvtk flye emu r-base
R
install.packages('ggplot2')
# you will need to pick a CRAN mirror here, I usually just pick one at random
install.packages('dplyr')
install.packages('BiocManager')
BiocManager::install('decontam')
```

## GENERAL TIPS
Some of these scripts will take a while to run, and will fail if they lose connection to the server, so I always run them like this:   

```bash
nohup sh script.sh input output > nohup_out &
```

Adding `nohup` before the command will tell it to ignore any hangup signals, like when the connection to the server drops  
Using `> nohup_out` will redirect all the stuff that would normally be printed to the screen into a file with that name  
The `&` at the end will tell it to run in the background so that you can do other things on the command line  

If you want to run multiple jobs at the same time, make sure you have unique names for the `names` and `nohup_out` files so that there is no chance of using the wrong file or overwriting the input/output mid-job. 

## DATA PROCESSING

### ONT FASTQ BACKUP
This is very specific to the CPG-IH data storage structure - it will take the data outputted by the onION, store the necessary files in our mediaflux backup, rename them if necessary, and generate a sharing link if requested.   

### ONT FASTQ MERGE AND RENAME
The `ont_combine_fastq.sh` script that will take the FASTQ files outputted by an ONT sequencer and do some processing steps that are sometimes necessary. The default data structure for FASTQ files outputted by ONT machines is a directory (eg. `fastq_pass`) containing a subdirectory for each sample named by the barcode detected (eg. `fastq_pass/barcode01/`). These subdirectories can contain a single `fastq.gz` file or multiple files.  
This script will detect if there are multiple files per sample and combine them to a single file per sample, renaming that file based on information provided.  
Four positional arguments are required:
1. The input directory containing the sample subdirectories
2. The output directory you want to create 
3. A two column TSV file listing, for each sample, the barcode name (eg. `barcode01`) and the new sample name to use (eg. `sampleA`). One sample per line. 
4. The output format desired - this can take one of two values:  
    * `file` - the output directory will contain a single `fastq.gz` file for each sample
    * `subdir` - the output directory will contain a subdirectory for each sample, each containing a single `fastq.gz` file
    * if unspecified, or an invalid foramt is specified, it will default to `file`

Example:
```bash
sh ~/Tools/CPGIH_Utility/ont_combine_fastq.sh inputdirectory outputdirectory renaming.tsv file
```

### TIDY PUBLIC HTML DATA
Another one specific to the CPG-IH data storage and sharing system, this time to tidy up the duplicated data that we use to generate the sharing links - deleting data that was shared more than 30 days ago.  
It takes one positional argument, the path to your `public_html/tmp` directory.  

Example:
```bash
bash ~/Tools/CPGIH_Utility/tidy_public_html.sh ~/public_html/tmp/
```

## DATA QUALITY CONTROL AND PRELIMINARY ANALYSIS

### ONT 16S AMPLICONS 
The `ont_qcemu.sh` script will take demultiplexed FASTQ reads and perform length filtering and taxonomic profiling.
Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the file before `.fastq.gz`), one per line
2. The directory where those FASTQ files are stored
3. The output directory you want to create  

If you want to change the length filtering cutoffs you can open the script in your text editor of choice and modify the `MINLEN` and `MAXLEN` variables at the top of the script. The defaults of these are 1400 and 1700 respectively.  

Example:
```bash
sh ~/Tools/CPGIH_Utility/ont_qcemu.sh names inputdirectory outputdirectory
```

The `barplots.R` script will make stacked barplots showing the 25 most abundant species identfied by Emu in the `ONT 16S DATA PROCESSING` section above.  
Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the file before `.fastq.gz`), one per line
2. The path to the Emu output generated in the previous step is stored (will likely be in `outputdirectory/EmuResults/`)
3. The output PDF you want to create  

Examples:
```bash
Rscript ~/Tools/CPGIH_Utility/barplots.R names emu-combined-abundance-species.tsv barplot.pdf

Rscript ~/Tools/CPGIH_Utility/barplots.R names emu-combined-abundance-species.tsv barplot.pdf 24 8
```

The first option will run the script at default, the second option will modify the width and height (in inches) of the output PDF he defaults are 12 and 8 respecively. If you want to modify the height or width, you will need to specify both - even if the other is the same as a default value.  

### ILLUMINA ISOLATE GENOMES
The `illumina_genomesqc.sh` script will take the input FASTQ reads, assign taxonomy to the reads using Kraken2, generate a draft assembly with shovill, and summarise the assembly quality (length, N50 etc.) and (depth of) coverage.  
Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the sample), one per line
2. The path to the directory containing the `fastq.gz` reads
3. The output directory you want to create 

Example:
```bash
bash ~/Tools/CPGIH_Utility/illumina_genomesqc.sh names inputdirectory outputdirectory
```

### ONT ISOLATE GENOMES
The `ont_genomesqc.sh` script works similar to the Illumina one. It will take the input FASTQ reads, assign taxonomy to the reads using Kraken2, generate a draft assembly with flye, and summarise the assembly quality (length, N50 etc.) and (depth of) coverage.  
Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the sample), one per line
2. The path to the directory containing the `fastq.gz` reads
3. The output directory you want to create 

Example:
```bash
bash ~/Tools/CPGIH_Utility/ont_genomesqc.sh names inputdirectory outputdirectory
```
:construction: TO DO LIST :construction:
### ILLUMINA METAGENOMES

### ONT METAGENOMES

### EMP V4 AMPLICONS