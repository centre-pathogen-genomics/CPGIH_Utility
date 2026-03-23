from setuptools import setup

setup(
    name="cpgih_utility",
    version="1.0.0",
    description="All the stuff needed for QC and prelim analysis of seq data generated in CPG-IH",
    author="Calum J. Walsh",

    scripts=[
        'Scripts/calcgenbankgenomesize.sh',
        'Scripts/decontam.R',
        'Scripts/illumina_genomesqc.sh',
        'Scripts/illumina_metagenomesqc.sh',
        'Scripts/illumina_seq',
        'Scripts/ont_backup_stats_share.sh',
        'Scripts/ont_combine_fastq.sh',
        'Scripts/ont_filter_fastq.py',
        'Scripts/ont_genomesqc.sh',
        'Scripts/ont_qcemu.sh',
        'Scripts/ont_seq',
        'Scripts/ont_seq_assembly',
        'Scripts/ont_seq_short',
        'Scripts/tidy_public_html.sh'
        ],
)