#!/usr/bin/env perl

# returns list of histidines based on input, structure, and propka
#
# 2016, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage: assignhis.pl -hsd reslist -hse reslist -hsp reslist -auto [pdbFile]\n";
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
my $hsd=undef;
my $hse=undef;
my $hsp=undef;
my $allhse=0;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-hsd") {
    shift @ARGV;
    $hsd=shift @ARGV;
  } elsif ($ARGV[0] eq "-hse") {
    shift @ARGV;
    $hse=shift @ARGV;
  } elsif ($ARGV[0] eq "-hsp") {
    shift @ARGV;
    $hsp=shift @ARGV;
  } elsif ($ARGV[0] eq "-allhse") {
    shift @ARGV;
    $allhse=1;
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

my @hislist=();
my %histype;

my $c=$mol->{chain}->[0];
foreach my $r ( @{$c->{res}} ) {
  my $rname=$r->{name};
  if ($rname eq "HIS" || $rname eq "HSD" || $rname eq "HSE" || $rname eq "HSP") {
    push(@hislist,$r->{num});
    if ($allhse) {
      $histype{$r->{num}}="HSE";
    } else {
      $histype{$r->{num}}=$rname;
    }
  } 
}

if (defined $hsd) {
  foreach my $f ( split(/[:=]/,$hsd) ) {
    if (exists $histype{$f}) {
      $histype{$f}="HSD";
    }
  }
}
if (defined $hse) {
  foreach my $f ( split(/[:=]/,$hse) ) {
    if (exists $histype{$f}) {
      $histype{$f}="HSE";
    }
  }
}
if (defined $hsp) {
  foreach my $f ( split(/[:=]/,$hsp) ) {
    if (exists $histype{$f}) {
      $histype{$f}="HSP";
    }
  }
}

foreach my $n ( @hislist ) {
  if ($histype{$n} eq "HIS") {
    $histype{$n}="HSE";
  }
} 
 
my @allhsd=();
my @allhse=();
my @allhsp=();

foreach my $n ( @hislist ) {
  if ($histype{$n} eq "HSD") {
    push(@allhsd,$n);
  } elsif ($histype{$n} eq "HSE") {
    push(@allhse,$n);
  } elsif ($histype{$n} eq "HSP") {
    push(@allhsp,$n);
  }
}

my @option=();
push(@option,sprintf "hsd=%s",join(":",@allhsd)) if ($#allhsd>=0);
push(@option,sprintf "hse=%s",join(":",@allhse)) if ($#allhse>=0);
push(@option,sprintf "hsp=%s",join(":",@allhsp)) if ($#allhsp>=0);

printf "%s\n",join(",",@option) if ($#option>=0);
