
# DELETE READ SHARING DIRECTORIES OLDER THAN 30 DAYS

TMPDIR=$1

find ${TMPDIR}/* -type d -ctime +30 | xargs rm -rf
