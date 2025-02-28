# CPG Innovation Hub Utility Scripts

This repository contains utility scripts for random stuff we do in CPG-IH.  
You can use the links immediately below here to more easily navigate to relevant sections

:construction: convert into wiki format

- [INSTALLATION AND SETUP](#installation-and-setup)
- [DOWNLOADING AND INSTALLING DATABASES](#downloading-and-installing-databases)
    - [KRAKEN2](#kraken2)
- [GENERAL TIPS](#general-tips)
- [DATA PROCESSING](#data-processing)
    - [ONT FASTQ BACKUP](#ont-fastq-backup)
    - [ONT FASTQ MERGE AND RENAME](#ont-fastq-merge-and-rename)
    - [TIDY PUBLIC HTML](#tidy-public-html)
- [DATA QUALITY CONTROL AND PRELIMINARY ANALYSIS](#data-quality-control-and-preliminary-analysis)
    - [ONT 16S AMPLICONS](#ont-16S-amplicons)
    - [ILLUMINA ISOLATE GENOMES](#illumina-isolate-genomes)
    - [ONT ISOLATE GENOMES](#ont-isolate-genomes)

## INSTALLATION AND SETUP
First you need to clone this repository to your local machine/server. You can do this to any directory, but for the purposes of this documentation we will do it in a new directory called `Tools/` in your home directory (`~`).  

```bash
# create directory to store repository locally
mkdir ~/Tools/

# clone repository to newly created directory
git clone https://github.com/centre-pathogen-genomics/CPGIH_Utility.git ~/Tools/CPGIH_Utility
```

If you are using the MDU servers is it easiest to just load my conda env - it should set all the paths correctly.  

```bash
conda activate /home/cwwalsh/miniconda3/envs/cpgih_utility
```

Alternatively you can make your own conda environment.  

:construction: add note here on installing conda locally  
:construction: remove decontam and barplots from pipelines as they are too difficult to reliably install in conda env - instead include the code as a tutorial to be performed manually  

Following these steps should give you one that does everything - you will need to install your own databases for [Kraken2](https://benlangmead.github.io/aws-indexes/k2) and [emu](https://github.com/treangenlab/emu) (see below for instructions) or modify the scripts to use existing ones.  

```bash
# create a conda environment into which you will install the software
conda create -n cpgih_utility -y

# activate the conda environment
conda activate cpgih_utility

# install the required software 
conda install -c bioconda kraken2 shovill seqkit csvtk flye emu -y
```

## DOWNLOADING AND INSTALLING DATABASES 

### KRAKEN2
The scripts for Quality Control of genomic and metagenomic data use Kraken2 to identify the microbial species present. For this you will need to install suitable databases.  
There are lot of Kraken2 databases to choose from [here](https://benlangmead.github.io/aws-indexes/k2) - generally those with a greater phylogenetic range (covering bacteria, archaea, fungi, protists etc.) will be larger and require more computing time and resources.  
For the majority of users interested in identify "non-weird" microbial isolates and profiling human microbiome data, the Standard-8 or PlusPF-8 will suffice.   
```bash
# make a database in your home directory to store the database
mkdir ~/kraken2_db

# download the database to this directory - this will take a while
curl -o ~/kraken2_db/k2_standard_08gb_20241228.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz

# extract the database
tar -xzf ~/kraken2_db/k2_standard_08gb_20241228.tar.gz -C ~/kraken2_db 

# remove the original download to save space
rm ~/kraken2_db/k2_standard_08gb_20241228.tar.gz 

# inpect the database to make sure everything is set up correctly
kraken2-inspect --db ~/kraken2_db | head

# you should get an output that looks something like this
# Database options: nucleotide db, k = 35, l = 31
# Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
# Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
# Total taxonomy nodes: 50914
# Table size: 1398394626
# Table capacity: 2000000000
# Min clear hash value = 16812150170552094720
100.00	1398394626	458617	R	1	root
 99.28	1388303726	356261	R1	131567	  cellular organisms
 91.48	1279314785	2721882	D	2	    Bacteria

# the final thing you will need to do is tell conda where to find this database
# you can do this by setting the KRAKEN_DEFAULT_DB variable 
# note that kraken2 doesn't like the '~' symbol in the path so we will fill that in
conda env config vars set KRAKEN2_DEFAULT_DB=$HOME/kraken2_db

# you should then see a message like this
To make your changes take effect please reactivate your environment

# so now we deactivate and reactivate the environment to make the changes
conda deactivate
conda activate cpgih_utility

# confirm the changes have taken effect
conda env config vars list

# you will get an output that looks something like this
# note the path will be different based on your system and username
KRAKEN2_DEFAULT_DB = /Users/cwwalsh/kraken2_db

# this will now be set every time you load the conda environment
```

:construction: fix QC scripts to use default kraken2 database on roosta

## GENERAL TIPS
Some of the scripts in this repository will take a while to run - if you are running these on a remote server without a job scheduler (eg. SLURM) then they will fail if you lose connection to the server.  
You can avoid that by using something like `tmux` or `screen`, but I always run them like this:   

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
The `ont_combine_fastq.sh` script will run some processing steps on the FASTQ files outputted by an ONT sequencer. The default data structure for FASTQ files outputted by ONT machines is a directory (eg. `fastq_pass`) containing a subdirectory for each sample named by the barcode detected (eg. `fastq_pass/barcode01/`). These subdirectories can contain a single `fastq.gz` file or multiple files.  
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
sh ~/Tools/CPGIH_Utility/Scripts/ont_combine_fastq.sh inputdirectory outputdirectory renaming.tsv file
```

### TIDY PUBLIC HTML
Another one specific to the CPG-IH data storage and sharing system, this time to tidy up the duplicated data that we use to generate the sharing links - deleting data that was shared more than 30 days ago.  
It takes one positional argument, the path to your `public_html/tmp` directory.  

Example:
```bash
bash ~/Tools/CPGIH_Utility/Scripts/tidy_public_html.sh ~/public_html/tmp/
```

## DATA QUALITY CONTROL AND PRELIMINARY ANALYSIS

### ONT 16S AMPLICONS 
The `ont_qcemu.sh` script will take demultiplexed FASTQ reads and perform length filtering and [taxonomic profiling](https://github.com/treangenlab/emu).  

Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the file before `.fastq.gz`), one per line
2. The directory where those FASTQ files are stored
3. The output directory you want to create  

If you want to change the length filtering cutoffs you can open the script in your text editor of choice and modify the `MINLEN` and `MAXLEN` variables at the top of the script. The defaults of these are 1400 and 1700 respectively.  

Example:
```bash
sh ~/Tools/CPGIH_Utility/Scripts/ont_qcemu.sh names inputdirectory outputdirectory
```

The `barplots.R` script will make stacked barplots showing the 25 most abundant species identfied in the `ONT 16S DATA PROCESSING` section above.  

Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the file before `.fastq.gz`), one per line
2. The path to the Emu output generated in the previous step is stored (will likely be in `outputdirectory/EmuResults/`)
3. The output PDF you want to create  

Examples:
```bash
Rscript ~/Tools/CPGIH_Utility/Scripts/barplots.R names emu-combined-abundance-species.tsv barplot.pdf

Rscript ~/Tools/CPGIH_Utility/Scripts/barplots.R names emu-combined-abundance-species.tsv barplot.pdf 24 8
```

The first option will run the script as default, the second option will modify the width and height (in inches) of the output PDF (the defaults are 12 and 8 respecively so this will double the width, useful if you have a lot of samples). If you want to modify the height or width, you will need to specify both - even if the other is the same as a default value.  

### ILLUMINA ISOLATE GENOMES
The `illumina_genomesqc.sh` script will take the input FASTQ reads, assign taxonomy to the reads using [Kraken2](https://github.com/DerrickWood/kraken2) and generate a draft assembly with [SPAdes](https://github.com/ablab/spades).  
Expected output is a directory with subdirectories for Kraken2 and Spades outputs - taxonomic profiles and assembled genomes respectively.  The user will likely be interested in the file `summary.tsv` which describes, for each sample: 
- the number of input FASTQ read pairs
- the distribution of the read lengths (min, average, max, N50) 
- the number of contigs in the final assembly
- total assembly length
- assembly N50
- top 3 most abundant species as detected by Kraken2 

When running the script, three positional arguments are required:
1. A file listing the basenames you want to include (the name of the sample), one per line
2. The path to the directory containing the `fastq.gz` reads
3. The output directory you want to create 

Example:
```bash
bash ~/Tools/CPGIH_Utility/Scripts/illumina_genomesqc.sh names inputdirectory outputdirectory
```

### ONT ISOLATE GENOMES
The `ont_genomesqc.sh` script works similar to the Illumina one. It will take the input FASTQ reads, assign taxonomy to the reads using [Kraken2](https://github.com/DerrickWood/kraken2) and generate a draft assembly with [flye](https://github.com/mikolmogorov/Flye).  
Expected output is, like the Illumina script above, a directory with subdirectories for Kraken2 and Flye outputs - taxonomic profiles and assembled genomes respectively and the `summary.tsv` file describing, for each sample: 
- the number of input FASTQ reads
- the distribution of the read lengths (min, average, max, N50) 
- the number of contigs in the final assembly
- total assembly length
- assembly N50
- top 3 most abundant species as detected by Kraken2 

Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the sample), one per line
2. The path to the directory containing the `fastq.gz` reads
3. The output directory you want to create 

Example:
```bash
bash ~/Tools/CPGIH_Utility/Scripts/ont_genomesqc.sh names inputdirectory outputdirectory
```
:construction: :construction: :construction: :construction: :construction:
### ILLUMINA METAGENOMES

### ONT METAGENOMES

### EMP V4 AMPLICONS