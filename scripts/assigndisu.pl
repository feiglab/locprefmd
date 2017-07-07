#!/usr/bin/env perl

# returns list of  disulfide bonds to be added
#
# 2016, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage: assigndisu.pl -disu reslist -auto [pdbFile]\n";
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
my $disu=undef;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-disu") {
    shift @ARGV;
    $disu=shift @ARGV;
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

my @cyslist=();
my %havecys;

my $c=$mol->{chain}->[0];
foreach my $r ( @{$c->{res}} ) {
  my $rname=$r->{name};
  if ($rname eq "CYS") {
    push(@cyslist,$r);
    $havecys{$r->{num}}=1;
  } 
}

foreach my $s ( @{$mol->{ssbond}} ) {
  $havecys{$s->{resnum1}}=2;
  $havecys{$s->{resnum2}}=2;
}

my @disulist=();
if (defined $disu) {
  foreach my $d ( split(/=/,$disu) ) {
    my @f=split(/:/,$d);
    if ($havecys{$f[0]} == 1 && $havecys{$f[1]} == 1) {
      $havecys{$f[0]}=2;
      $havecys{$f[1]}=2;
      push(@disulist,$d);
    }
  }
}

if ($auto) {
  $mol->findSSBonds();
  foreach my $s ( @{$mol->{ssbond}}) {
    if ($havecys{$s->{resnum1}} == 1 && $havecys{$s->{resnum2}} == 1) {
      $havecys{$s->{resnum1}}=2;
      $havecys{$s->{resnum2}}=2;
      push(@disulist,sprintf("%d:%d",$s->{resnum1},$s->{resnum2}));
    }
  }
}

printf "-ssbond %s\n",join("=",@disulist) if ($#disulist>=0);
