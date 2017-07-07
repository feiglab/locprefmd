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
my $analyze;

my $cmol=Molecule::new();
$cmol->readPDB($pdb);
&Analyze::phipsi($cmol);

my @needtofixcis=();
my $c=$cmol->{chain}->[0];
for (my $ir=0; $ir<$#{$c->{res}}; $ir++) {
  my $r=$c->{res}->[$ir];
  if ($r->{omega}>-90 && $r->{omega}<90) {
    if ($r->{name} ne "PRO" && $c->{res}->[$ir+1]->{name} ne "PRO") {
      push(@needtofixcis,$ir);
    }
  }
}

my $sicho=SICHO::new(gridsize => 100, offsetx  => 50.0, offsety  => 50.0,
                     offsetz  => 50.0, resolution => 1.0, intflag    => 0);
$sicho->genSimpleFromAllAtom($cmol,ca=>1);

my $seq=Sequence::new($cmol);

my $rebmol=Molecule::new();
$rebmol->rebuildFromSICHO($seq,$sicho,undef,undef,1);
$rebmol->resetValidResidues(0,1);

$c=$rebmol->{chain}->[0];
my $atom=$c->{atom};
my $res=$c->{res};
foreach my $t ( @needtofixcis ) {
  my $r=$res->[$t];
#  printf STDERR "%s::%s::%s inx %d\n",$r->{chain},$r->{num},$r->{name},$t;
  $r->{valid}=1;
  for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
    $atom->[$ia]->{valid}=1;
  }
}

my $smol=Molecule::new();
$smol=$rebmol->clone(1);
$cmol->merge($smol);

$cmol->writePDB("-","GENERIC");




