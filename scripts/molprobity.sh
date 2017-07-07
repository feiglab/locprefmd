#!/bin/bash
source $MOLPROBITY/build/setpaths.sh
convpdb.pl -out generic $1 | pdbelement.pl > /tmp/$$.pdb
$MOLPROBITY/cmdline/molprobity /tmp/$$.pdb | tail -1 | awk -F: '{print $32}'
/bin/rm -f /tmp/$$.pdb
