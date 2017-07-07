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
  if (($r->{omega}>-90 && $r->{omega}<90) || 1) {
    if ($r->{name} ne "PRO" && $c->{res}->[$ir+1]->{name} ne "PRO") {
      push(@needtofixcis,$r);
    }
  }
}

foreach my $r ( @needtofixcis ) {
  printf "CONS DIHE %s %d CA %s %d C %s %d N %s %d CA FORCE 10.0 MIN 180 WIDTH 10 PERIOD 0\n",
     $r->{seg},$r->{num},$r->{seg},$r->{num},$r->{seg},$r->{num}+1,$r->{seg},$r->{num}+1;
}
