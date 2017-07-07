#!/usr/bin/env perl

# returns list of patches to protonated ASP/GLU based on input, structure, and propka
#
# 2016, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage: assignaspglu.pl -asp reslist -glu reslist -auto [pdbFile]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Molecule;

my $fname;

my $auto=0;
my $asp=undef;
my $glu=undef;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-asp") {
    shift @ARGV;
    $asp=shift @ARGV;
  } elsif ($ARGV[0] eq "-glu") {
    shift @ARGV;
    $glu=shift @ARGV;
  } elsif ($ARGV[0] eq "-auto") {
    shift @ARGV;
    $auto=1;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $fname = shift @ARGV;
    $done=1;
  }
}

my $mol=Molecule::new();
$mol->readPDB($fname);

my @aspglulist=();
my %aspglutype;

my $c=$mol->{chain}->[0];
foreach my $r ( @{$c->{res}} ) {
  my $rname=$r->{name};
  if ($rname eq "ASP" || $rname eq "GLU") {
    push(@aspglulist,$r->{num});
    $aspglutype{$r->{num}}=$rname;
  } 
}

my @patches=();
if (defined $asp) {
  foreach my $f ( split(/[:=]/,$asp) ) {
    if (exists $aspglutype{$f} && $aspglutype{$f} eq "ASP") {
      push(@patches,sprintf("ASPP:%d",$f));
    }
  }
}
if (defined $glu) {
  foreach my $f ( split(/[:=]/,$glu) ) {
    if (exists $aspglutype{$f} && $aspglutype{$f} eq "GLU") {
      push(@patches,sprintf("GLUP:%d",$f));
    }
  }
}

printf "patch=%s\n",join("_",@patches) if ($#patches>=0);
