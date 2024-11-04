# Centre for Pathogen Genomics: Innovation Hub Utility Scripts

This repository contains utility scripts for random stuff we do in CPG-IH.  

### INSTALLATION AND SETUP
First you need to clone this repository t oyour local machine/server. You can do this anywhere, but for the purposes for this documentation we will do it in a new directory called `Tools/` in your home directory (`~`) to make it easy.  

```bash
mkdir ~/Tools/
cd ~/Tools/
git clone https://github.com/centre-pathogen-genomics/CPGIH_Utility.git 
```

If you are using the MDU servers is it easiest to just load my conda env - should set all the paths correctly.  

```bash
conda activate /home/cwwalsh/miniconda3/envs/cpgih_utility
```

Alternatively you can make your own conda environment. Following these steps should give you on that does everything - you will need to install your own databases for `kraken2` and `emu` though.  

```bash
conda create -n cpgih_utility -y
conda activate cpgih_utility
conda install -c bioconda kraken2 shovill seqkit csvtk flye emu r-base
R
install.packages('tidyverse')
# you will need to pick a CRAN mirror here, I usually just pick Ireland 

```

### GENRAL TIPS
Some of these scripts will take a while to run (I usually let these run overnight), and will fail if they loses connection to the server, so I always run them like this if I am using a system without a job manager like SLURM:  

```bash
nohup sh script.sh input output > nohup_out &
```

Adding `nohup` before the command will tell it to ignore any hangup signals, like when the connection to the server drops  
`> nohup_out` will tell it to redirect all the stuff that would normally be printed to the screen into a file with that name  
the `&` at the end will tell it to run in the background so that you can do other things on the command line  

If you want to run multiple jobs at the same time, make sure you have unique names for the `names` and `nohup_out` filesso that there is no chance of using the wrong file or overwriting the input/output mid-job. 

### ONT FASTQ BACKUP
This is very specific to the CPG-IH data storage structure - it will take the data outputted by the minION(s), store the necessary files in our mediaflux backup, rename them if necessary, and generate a sharing link if requested.   

### ONT FASTQ MERGE AND RENAME



### ONT 16S DATA PROCESSING
The `ont_qcemu.sh` script will take demultiplexed FASTQ reads and perform length filtering and taxonomic profiling.
Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the file before `.fastq.gz`), one per line
2. The directory where those FASTQ files are stored
3. The output directory you want to create  

If you want to change the length filtering cutoffs you can open the script in your text editor of choice and modify the `MINLEN` and `MAXLEN` variables at the top of the script. The defaults of these are 1400 and 1700 respectively.  

Example:
```bash
sh ~/Tools/CPGIH_Utility/ont_qcemu.sh names inputdirectory outputdirectory > nohup_out &
```

### ONT 16S DATA PLOTTING
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

### ILLUMINA ISOLATE GENOME DATA PROCESSING
The `illumina_genomesqc.sh` script will take the input FASTQ reads, assign taxonomy to the reads using Kraken2, generate a draft assembly with shovill, and summarise the assembly quality (length, N50 etc.) and (depth of) coverage.  
Three positional arguments are required:
1. A file listing the basenames you want to include (the name of the sample), one per line
2. The path to the directory containing the `fastqgz` reads
3. The output directory you want to create 

Example:
```bash
~/Tools/CPGIH_Utility/ont_genomesqc.sh names inputdirectory outputdirectory
```

### ONT ISOLATE GENOME DATA PROCESSING
### TIDY PUBLIC HTML DATA

