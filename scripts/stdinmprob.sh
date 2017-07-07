#!/bin/bash

rnd=$RANDOM
host=$HOST
tmpname="$host-$rnd-$$"
convpdb.pl -out generic | pdbelement.pl > /tmp/$tmpname.pdb
source $MOLPROBITY/build/setpaths.sh
$MOLPROBITY/cmdline/molprobity /tmp/$tmpname.pdb | tail -1 | awk -F: '{print $32}' 2>/dev/null
rm -f /tmp/$tmpname.pdb
