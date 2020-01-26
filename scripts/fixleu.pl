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

my $pdb=shift @ARGV;

my $cmol=Molecule::new();
$cmol->readPDB($pdb,ignoreseg=>0);

my @reslist=();
my $c=$cmol->{chain}->[0];
for (my $ir=0; $ir<$#{$c->{res}}; $ir++) {
  my $r=$c->{res}->[$ir];
  if ($r->{name} eq "LEU") {
    my %have;
    my %haveinx;
    my $atom=$c->{atom};
    for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
      $have{$atom->[$ia]->{atomname}}=$atom->[$ia];
    }
    if (defined $have{CA} && defined $have{CB} && 
        defined $have{CG} && defined $have{CD1} && 
        defined $have{CD2}) {
       my $dihed1=&GenUtil::dihedral($have{CA},$have{CB},$have{CG},$have{CD1});
       my $dihed2=&GenUtil::dihedral($have{CA},$have{CB},$have{CG},$have{CD2});
       my $diff=$dihed2-$dihed1;
       $diff+=360 if ($diff<0);
       if ($diff>180 && $diff<300) {
         printf STDERR "%s fixing LEU %d\n",$pdb,$r->{num};
         my $tmpx=$have{CD1}->{xcoor};
         my $tmpy=$have{CD1}->{ycoor};
         my $tmpz=$have{CD1}->{zcoor};
         $have{CD1}->{xcoor}=$have{CD2}->{xcoor};
         $have{CD1}->{ycoor}=$have{CD2}->{ycoor};
         $have{CD1}->{zcoor}=$have{CD2}->{zcoor};
         $have{CD2}->{xcoor}=$tmpx;
         $have{CD2}->{ycoor}=$tmpy;
         $have{CD2}->{zcoor}=$tmpz;

         if (exists $have{HD11} && exists $have{HD21}) {
           my $tmpx=$have{HD11}->{xcoor};
           my $tmpy=$have{HD11}->{ycoor};
           my $tmpz=$have{HD11}->{zcoor};
           $have{HD11}->{xcoor}=$have{HD21}->{xcoor};
           $have{HD11}->{ycoor}=$have{HD21}->{ycoor};
           $have{HD11}->{zcoor}=$have{HD21}->{zcoor};
           $have{HD21}->{xcoor}=$tmpx;
           $have{HD21}->{ycoor}=$tmpy;
           $have{HD21}->{zcoor}=$tmpz;
         }
         if (exists $have{HD12} && exists $have{HD22}) {
           my $tmpx=$have{HD12}->{xcoor};
           my $tmpy=$have{HD12}->{ycoor};
           my $tmpz=$have{HD12}->{zcoor};
           $have{HD12}->{xcoor}=$have{HD22}->{xcoor};
           $have{HD12}->{ycoor}=$have{HD22}->{ycoor};
           $have{HD12}->{zcoor}=$have{HD22}->{zcoor};
           $have{HD22}->{xcoor}=$tmpx;
           $have{HD22}->{ycoor}=$tmpy;
           $have{HD22}->{zcoor}=$tmpz;
         }
         if (exists $have{HD13} && exists $have{HD23}) {
           my $tmpx=$have{HD13}->{xcoor};
           my $tmpy=$have{HD13}->{ycoor};
           my $tmpz=$have{HD13}->{zcoor};
           $have{HD13}->{xcoor}=$have{HD23}->{xcoor};
           $have{HD13}->{ycoor}=$have{HD23}->{ycoor};
           $have{HD13}->{zcoor}=$have{HD23}->{zcoor};
           $have{HD23}->{xcoor}=$tmpx;
           $have{HD23}->{ycoor}=$tmpy;
           $have{HD23}->{zcoor}=$tmpz;
         }
       }
    }
  }
}

$cmol->writePDB("-","GENERIC");

