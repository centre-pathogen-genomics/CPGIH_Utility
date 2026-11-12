#!/usr/bin/env bash

set -eu -o pipefail

# VERSION: one repo-wide number, stored in version.txt next to this script.
# Bump it (MAJOR = breaking change to inputs/outputs, MINOR = new feature, PATCH = bug fix) whenever any script changes
if [ "$1" = "--version" ] || [ "$1" = "-v" ]
then
    VERSIONFILE="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/version.txt"
    echo "$(basename "$0") $(cat "${VERSIONFILE}" 2>/dev/null || echo 'unknown (version.txt not found)')"
    exit 0
fi

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

sed '1d' assembly_summary_genbank.txt | \
    cut -f 7,26 assembly_summary_genbank.txt | \
    csvtk summary -t -f genome_size:mean -f genome_size:median -f genome_size:stdev -g species_taxid | \
    taxonkit lineage -i 1 | \
    taxonkit reformat -i 5 | \
    csvtk cut -t -f 1,2,3,4,6 > taxid_genomelen_lineage.txt

rm -f assembly_summary_genbank.txt
