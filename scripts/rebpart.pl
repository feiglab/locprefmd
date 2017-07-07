#!/usr/bin/env perl

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Molecule;
use Analyze;
use Sequence;
use SICHO;

my $pdb=shift @ARGV;
my $inx=shift @ARGV;
my $fixca=shift @ARGV;

$fixca=0 if (!defined $fixca);

my @finx=split(/[,:=]/,$inx);

my $mol=Molecule::new();
$mol->readPDB($pdb);

my @needtofix=();
my $c=$mol->{chain}->[0];
for (my $ir=0; $ir<$#{$c->{res}}; $ir++) {
  my $r=$c->{res}->[$ir];
  my $found=0;
  foreach my $f ( @finx ) {
    $found=1 if ($r->{num}>=$f-1 && $r->{num}<=$f+1);
  } 
  if ($found) {
    push(@needtofix,$ir);
  }
}

my $sicho=SICHO::new(gridsize => 100, offsetx  => 50.0, offsety  => 50.0,
                     offsetz  => 50.0, resolution => 1.0, intflag    => 0);
$sicho->genSimpleFromAllAtom($mol,ca=>$fixca);

my $seq=Sequence::new($mol);

my $rebmol=Molecule::new();
$rebmol->rebuildFromSICHO($seq,$sicho,undef,undef,$fixca);
$rebmol->resetValidResidues(0,1);

$c=$rebmol->{chain}->[0];
my $atom=$c->{atom};
my $res=$c->{res};
foreach my $t ( @needtofix ) {
  my $r=$res->[$t];
#  printf STDERR "%s::%s::%s inx %d\n",$r->{chain},$r->{num},$r->{name},$t;
  $r->{valid}=1;
  for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
    $atom->[$ia]->{valid}=1;
  }
}

my $smol=Molecule::new();
$smol=$rebmol->clone(1);
$mol->merge($smol);

$mol->writePDB("-","GENERIC");





