#!/bin/bash
# carry out local structure refinement for a single chain
# Michael Feig, Michigan State University, 2016

cpus=8
fixca=1
minsteps=200
mdsteps=5000
savefrq=1000
tmpdir=/tmp
hsdlist=
hselist=
hsplist=
aspplist=
gluplist=
disulist=
hisauto=1
disuauto=1

export PATH=$PATH:$LOCPREFMD/scripts

while :; do
  case $1 in 
    -c|--cpus) 
      if [ -n "$2" ]; then
        cpus=$2
        shift
      fi
      ;;
    -ca|--fixca)
      fixca=1
      ;;
    --dontfixca)
      fixca=0
      ;;
    --mdsteps)
      if [ -n "$2" ]; then
        mdsteps=$2
        shift
      fi
      ;;
    --minsteps)
      if [ -n "$2" ]; then
        minsteps=$2
        shift
      fi
      ;;
    --mdsavefrq)
      if [ -n "$2" ]; then
        savefrq=$2
        shift
      fi
      ;;
    --tmpdir)
      if [ -n "$2" ]; then
        tmpdir=$2
        shift
      fi
      ;;
    --disu) 
      if [ -n "$2" ]; then
        disulist=$2
        shift
      fi
      ;;
    --aspp) 
      if [ -n "$2" ]; then
        aspplist=$2
        shift
      fi
      ;;
    --glup) 
      if [ -n "$2" ]; then
        gluplist=$2
        shift
      fi
      ;;
    --hsp) 
      if [ -n "$2" ]; then
        hsplist=$2
        shift
      fi
      ;;
    --hse) 
      if [ -n "$2" ]; then
        hselist=$2
        shift
      fi
      ;;
    --hsd) 
      if [ -n "$2" ]; then
        hsdlist=$2a
        shift
      fi
      ;;
    --hisauto) 
      hisauto=1
      ;;
    --disuauto)
      disuauto=1
      ;;
    --)
      shift
      break
      ;;
    *)
      break
  esac
  shift
done
          
inp=$1

hostname=`uname -n`
ttag=$hostname-$$
tag="${ttag,,}"

pwd=`pwd`

/bin/rm -rf $tmpdir/$tag
mkdir -p $tmpdir/$tag
cd $tmpdir/$tag

topfile=$LOCPREFMD/data/xtra_top_all36_prot.rtf
parfile=$LOCPREFMD/data/xtrabonds2cmap_par_all36_prot.prm
ffpar="dielec=rdie,epsilon=4,param=x,xpar=$parfile,xtop=$topfile"

convpdb.pl -cleanaux -nsel protein $pwd/$inp > $tag.protein.pdb

hsdopt=
hseopt=
hspopt=
hisautoopt=

if [ $hsdlist ]; then
  hsdopt="-hsd $hsdlist"
fi
if [ $hselist ]; then
  hseopt="-hse $hselist"
fi
if [ $hsplist ]; then
  hspopt="-hsp $hsplist"
fi
if [ $hisauto -eq 1 ]; then
  hisautoopt="-auto"
fi

hisoptions=`assignhis.pl $hsdopt $hseopt $hspopt $hisautoopt $tag.protein.pdb`

if [ $hisoptions ]; then
  ffpar="$ffpar,$hisoptions"
fi

aspopt=
gluopt=

if [ $aspplist ]; then
  aspopt="-asp $aspplist";
fi
if [ $glupplist ]; then
  gluopt="-glu $gluplist";
fi

aspgluoptions=`assignaspglu.pl $aspopt $gluopt $tag.protein.pdb`

if [ $aspgluoptions ]; then
  ffpar="$ffpar,$aspgluoptions"
fi

ssbond=
ssauto=
if [ $disulist ]; then
  ssbond="-disu $disulist"
fi
if [ $disuauto -eq 1 ]; then
  ssauto="-auto"
fi

disuoptions=`assigndisu.pl $ssbond $ssauto $tag.protein.pdb`

chain=`head -1 $tag.protein.pdb | cut -c22`;

if [ "$chain" != " " ]; then
  convpdb.pl -nsel ${chain}: $tag.protein.pdb | fixleu.pl | convpdb.pl -setchain "X" | convpdb.pl $disuoption -segnames > $tag.init.pdb
else
  convpdb.pl -setchain "X" $tag.protein.pdb | fixleu.pl | convpdb.pl $disuoptions -segnames > $tag.init.pdb
fi

molprobity.sh $tag.init.pdb > $tag.init.mpscore &

for seg in `convpdb.pl -listseg $tag.init.pdb`; do
  convpdb.pl -nsel ${seg}: -readseg $tag.init.pdb | complete.pl | convpdb.pl -setseg $seg -setchain " " > $tag.$seg.pdb
  convpdb.pl -out generic $tag.$seg.pdb > $tag.$seg.clean.pdb
  convpdb.pl -out generic $tag.$seg.clean.pdb | convpdb.pl -out generic -setseg $seg > $tag.$seg.zero.pdb
  fixcis.pl $tag.$seg.clean.pdb | convpdb.pl -out generic -setseg $seg > $tag.$seg.fix.pdb
  fixciscons.pl $tag.$seg.zero.pdb >> $tag.consdihe.custom
  /bin/cat $tag.$seg.zero.pdb | egrep -v "(TER|END)" >> $tag.clean.pdb
  /bin/cat $tag.$seg.fix.pdb | egrep -v "(TER|END)" >> $tag.fix.pdb
done

molprobity.sh $tag.clean.pdb > $tag.clean.mpscore &

minCHARMM.pl -nochain -par sdsteps=0,minsteps=0,$ffpar -log $tag.m0log $tag.clean.pdb > $tag.clean.cinp.pdb
minCHARMM.pl -nochain -par sdsteps=0,minsteps=0,$ffpar -log $tag.fix.m0log $tag.fix.pdb > $tag.fix.cinp.pdb

minCHARMM.pl -nochain -custom $tag.consdihe.custom -cons ca $tag.clean.cinp.pdb 0:999_0.1 -par minsteps=$minsteps,$ffpar -log $tag.m1log $tag.fix.cinp.pdb > $tag.min1.pdb

xtag=1
inp=$tag.min1.pdb
caforce=0.1

checkring.pl $inp > $tag.checkring.out 2>&1
reslist=`cat $tag.checkring.out | awk '/CROSS/ {print $3}' | sort -u | xargs echo | sed -e "s/ /:/g"`
/bin/rm -f $tag.checkring.out

while [ "$reslist" -a $xtag -lt 6 ]; do
  convpdb.pl -setall -setchain " " -out generic -setseg PRO0 $inp | grep -v TER > $tag.tinp.pdb
  rebpart.pl $tag.tinp.pdb $reslist | convpdb.pl -setchain " " -setseg "PROX" | minCHARMM.pl -nochain -custom $tag.consdihe.custom -cons ca $tag.clean.cinp.pdb 0:999_$caforce -par minsteps=$minsteps,$ffpar -log "$tag.m1x${xtag}.log" > $tag.min1x$xtag.pdb
  /bin/rm -f $tag.tinp.pdb

  inp="$tag.min1x$xtag.pdb"

  checkring.pl $inp > $tag.checkring.out 2>&1
  reslist=`cat $tag.checkring.out | awk '/CROSS/ {print $3}' | sort -u | xargs echo | sed -e "s/ /:/g"`
  /bin/rm -f $tag.checkring.out

  xtag=$[$xtag+1]
  caforce=`echo $caforce | awk '{print $1/2.0}'`
done

/bin/rm -f $tag.min1y.pdb
/bin/ln -s $inp $tag.min1y.pdb

minCHARMM.pl -nochain -custom $tag.consdihe.custom -cons ca $tag.clean.cinp.pdb 0:999_0.2 -par minsteps=$minsteps,$ffpar -log $tag.m2log $tag.min1y.pdb > $tag.min2.pdb

minCHARMM.pl -nochain -custom $tag.consdihe.custom -cons ca $tag.clean.cinp.pdb 0:999_0.5 -par minsteps=$minsteps,$ffpar -log $tag.m3log $tag.min2.pdb > $tag.min3.pdb

inp=$tag.min3.pdb

reslist=`mprobrotamers.sh $inp | grep OUTLIER | awk '{print $2}' | xargs echo | sed -e "s/ /:/g"`

if [ "$reslist" ]; then
  convpdb.pl -setall -setchain " " -out generic -setseg PRO0 $inp | grep -v TER > $tag.tinp.pdb
  rebpart.pl $tag.tinp.pdb $reslist 1 | convpdb.pl -setchain X | convpdb.pl -segnames | minCHARMM.pl -nochain -custom $tag.consdihe.custom -cons ca $tag.clean.cinp.pdb 0:999_0.5 -par minsteps=$minsteps,$ffpar -log $tag.m3xlog > $tag.min3x.pdb
  /bin/rm -f $tag.tinp.pdb
else
  /bin/rm -f $tag.min3x.pdb
  /bin/ln -s $inp $tag.min3x.pdb
fi

minCHARMM.pl -nochain -custom $tag.consdihe.custom -cons ca $tag.clean.cinp.pdb 0:999_1 -par minsteps=$minsteps,$ffpar -log $tag.m4log $tag.min3x.pdb > $tag.min4.pdb

mkdir $tag.md{1,2,3,4}
cd $tag.md1
mdCHARMM.pl -nochain -custom ../$tag.consdihe.custom -cons ca ../$tag.clean.cinp.pdb 0:999_1 -par dynsteps=$mdsteps,dyntemp=20,dynber,dynbertc=0.1,dynoutfrq=$savefrq,dyntstep=0.001,$ffpar -log $tag.md1log -final ../$tag.md1.pdb -trajout ../$tag.md1.dcd ../$tag.min4.pdb &

cd ../$tag.md2
mdCHARMM.pl -nochain -custom ../$tag.consdihe.custom -cons ca ../$tag.clean.cinp.pdb 0:999_1 -par dynsteps=$mdsteps,dyntemp=100,dynber,dynbertc=0.1,dynoutfrq=$savefrq,dyntstep=0.001,$ffpar -log $tag.md2log -final ../$tag.md2.pdb -trajout ../$tag.md2.dcd ../$tag.min4.pdb &

cd ../$tag.md3
mdCHARMM.pl -nochain -custom ../$tag.consdihe.custom -cons ca ../$tag.clean.cinp.pdb 0:999_0.5 -par dynsteps=$mdsteps,dyntemp=150,dynber,dynbertc=0.1,dynoutfrq=$savefrq,dyntstep=0.001,$ffpar -log $tag.md3log -final ../$tag.md3.pdb -trajout ../$tag.md3.dcd ../$tag.min4.pdb &

cd ../$tag.md4
if [ $fixca -eq 0 ]; then
  mdCHARMM.pl -nochain -custom ../$tag.consdihe.custom -par dynsteps=$mdsteps,dyntemp=150,dynber,dynbertc=0.1,dynoutfrq=$savefrq,dyntstep=0.001,$ffpar -log $tag.md4log -final ../$tag.md4.pdb -trajout ../$tag.md4.dcd ../$tag.min4.pdb &
else
  mdCHARMM.pl -nochain -custom ../$tag.consdihe.custom -cons ca ../$tag.clean.cinp.pdb 0:999_2 -par dynsteps=$mdsteps,dyntemp=150,dynber,dynbertc=0.1,dynoutfrq=$savefrq,dyntstep=0.001,$ffpar -log ../$tag.md4log -final ../$tag.md4.pdb -trajout ../$tag.md4.dcd ../$tag.min4.pdb &
fi

cd ..
wait
imp=`cat $tag.init.mpscore`
impc=`cat $tag.clean.mpscore`

/bin/rm -rf $tag.optmd
for n in 1 2 3 4; do
  processDCD.pl -ensdir $tag.optmd -ens md $tag.min1.pdb $tag.md$n.dcd
done

if [ $fixca -eq 0 ]; then
  ensrun.pl -nocompress -dir $tag.optmd -cpus $cpus -new min md minCHARMM.pl -nochain -custom `pwd`/$tag.consdihe.custom -cons ca self 0:999_0.01 -par minsteps=$minsteps,$ffpar 
  else
  ensrun.pl -nocompress -dir $tag.optmd -cpus $cpus -new min md minCHARMM.pl -nochain -custom `pwd`/$tag.consdihe.custom -cons ca self 0:999_1 -par minsteps=$minsteps,$ffpar 
fi

checkin.pl -dir $tag.optmd min $tag.min*.pdb
ensrun.pl -dir $tag.optmd -cpus $cpus -set molprobity:1 min stdinmprob.sh > /dev/null 2>&1
ensrun.pl -dir $tag.optmd -cpus $cpus -set irmsdca:1 min rms.pl -fit -nowarn -out CA `pwd`/$tag.clean.cinp.pdb
ensrun.pl -dir $tag.optmd -cpus $cpus -set irmsdall:1 min rms.pl -fit -nowarn -out all `pwd`/$tag.clean.cinp.pdb
ensrun.pl -dir $tag.optmd -set molprobity:1 min stdinmprob.sh  > /dev/null 2>&1
ensrun.pl -dir $tag.optmd -set molprobity:1 min stdinmprob.sh  > /dev/null 2>&1

if [ $fixca -eq 0 ]; then
  best=`ensfiles.pl -prop molprobity,irmsdca,irmsdall -dir $tag.optmd min | awk '{print $1,$2+$3*0.001,$2,$3,$4}' | sort -k 2 -n | head -1` 
else 
  best=`ensfiles.pl -prop molprobity,irmsdca,irmsdall -dir $tag.optmd min | awk '{print $1,$2+$3*0.1,$2,$3,$4}' | sort -k 2 -n | head -1` 
fi
read -a results <<< $best
fmp=${results[2]}
irmsca=${results[3]}
irmssc=${results[4]}

convpdb.pl -setall -setchain "$chain" ${results[0]} -out generic > $tag.final.pdb

echo "REMARK initial molprobity $imp"
echo "REMARK initial/complete molprobity $impc"
echo "REMARK final molprobity $fmp"
echo "REMARK CA-RMSD from-initial $irmsca"
echo "REMARK SC-RMSD from-initial $irmssc"

if [ `grep ATOM $tag.final.pdb | wc -l` -lt 1 ]; then
  echo "ERROR: refinement unsuccessful"
else 
  cat $tag.final.pdb
fi

/bin/rm -rf $tmpdir/$tag
