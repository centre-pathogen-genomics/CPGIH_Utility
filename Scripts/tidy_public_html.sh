#!/bin/bash

# VERSION: one repo-wide number, stored in version.txt next to this script.
# Bump it (MAJOR = breaking change to inputs/outputs, MINOR = new feature, PATCH = bug fix) whenever any script changes
if [ "$1" = "--version" ] || [ "$1" = "-v" ]
then
    VERSIONFILE="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/version.txt"
    echo "$(basename "$0") $(cat "${VERSIONFILE}" 2>/dev/null || echo 'unknown (version.txt not found)')"
    exit 0
fi

# DELETE READ SHARING DIRECTORIES OLDER THAN 30 DAYS

TMPDIR=$1

find ${TMPDIR}/* -type d -ctime +30 | xargs rm -rf
