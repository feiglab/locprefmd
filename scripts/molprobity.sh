#!/bin/bash

source $MOLPROBITY/build/setpaths.sh

TMPDIR=/tmp/molp.$$
cwd=$(pwd)
mkdir -p $TMPDIR 

convpdb.pl -out generic $1 | pdbelement.pl > $TMPDIR/$$.pdb
cd $TMPDIR

$MOLP/cmdline/oneline-analysis . | tail -1 | awk -F: '{print $32}'

cd $cwd
/bin/rm -rf $TMPDIR
