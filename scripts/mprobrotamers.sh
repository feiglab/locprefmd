#!/bin/bash
source $MOLPROBITY/build/setpaths.sh
convpdb.pl -out generic $1 | pdbelement.pl > /tmp/$$.pdb
echo "Rotamer outliers"
$MOLPROBITY/build/bin/phenix.rotalyze data_version=8000 /tmp/$$.pdb | grep OUTLIER
rm -f /tmp/$$.pdb
