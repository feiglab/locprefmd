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

my $fname;

$fname = shift @ARGV;

my $mol=Molecule::new();
$mol->readPDB($fname,readseg=>1);


foreach my $c ( @{$mol->{chain}} ) {
  my $a=$c->{atom};
  foreach my $r ( @{$c->{res}}) {
    my $cx=0.0;
    my $cy=0.0;
    my $cz=0.0;
    my $nn=0;
    my $czx=0.0;
    my $czy=0.0;
    my $czz=0.0;
    my $cdx=0.0;
    my $cdy=0.0;
    my $cdz=0.0;

    if ($r->{name} eq "PHE" || $r->{name} eq "TYR") {
      for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
        my $n=$a->[$i]->{atomname};
        if ($n eq "CD1" || $n eq "CD2" || $n eq "CE1" || $n eq "CE2" || $n eq "CZ") {
           $cx+=$a->[$i]->{xcoor};
           $cy+=$a->[$i]->{ycoor};
           $cz+=$a->[$i]->{zcoor};
           $nn++;
           if ($n eq "CZ") {
             $czx=$a->[$i]->{xcoor};
             $czy=$a->[$i]->{ycoor};
             $czz=$a->[$i]->{zcoor};
           }
           if ($n eq "CD1") {
             $cdx=$a->[$i]->{xcoor};
             $cdy=$a->[$i]->{ycoor};
             $cdz=$a->[$i]->{zcoor};
           }
        } 
      }
      $cx/=$nn;
      $cy/=$nn;
      $cz/=$nn;

      my $dcx=$cx-$czx;
      my $dcy=$cy-$czy;
      my $dcz=$cz-$czz;
      my $dc=sqrt($dcx*$dcx+$dcy*$dcy+$dcz*$dcz);
      $dcx/=$dc;
      $dcy/=$dc;
      $dcz/=$dc;

      my $ddx=$cx-$cdx;
      my $ddy=$cy-$cdy;
      my $ddz=$cz-$cdz;
      my $dd=sqrt($ddx*$ddx+$ddy*$ddy+$ddz*$ddz);
      $ddx/=$dd;
      $ddy/=$dd;
      $ddz/=$dd;

      my $xx=$dcy*$ddz-$dcz*$ddy;
      my $xy=-$dcx*$ddz+$dcz*$ddx;
      my $xz=$dcx*$ddy-$dcy*$ddx;
      my $d=sqrt($xx*$xx+$xy*$xy+$xz*$xz);
      $xx/=$d;
      $xy/=$d;
      $xz/=$d;
 
      &checkIntersect($mol,$r,$cx,$cy,$cz,$xx,$xy,$xz);
    } elsif ($r->{name} eq "TRP") {
      for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
        my $n=$a->[$i]->{atomname};
        if ($n eq "CE2" || $n eq "CD2" || $n eq "CE3" || $n eq "CZ3" || $n eq "CZ2" || $n eq "CH2") {
           $cx+=$a->[$i]->{xcoor};
           $cy+=$a->[$i]->{ycoor};
           $cz+=$a->[$i]->{zcoor};
           $nn++;
           if ($n eq "CE3") {
             $czx=$a->[$i]->{xcoor};
             $czy=$a->[$i]->{ycoor};
             $czz=$a->[$i]->{zcoor};
           }
           if ($n eq "CH2") {
             $cdx=$a->[$i]->{xcoor};
             $cdy=$a->[$i]->{ycoor};
             $cdz=$a->[$i]->{zcoor};
           }
        } 
      }
      $cx/=$nn;
      $cy/=$nn;
      $cz/=$nn;

      my $dcx=$cx-$czx;
      my $dcy=$cy-$czy;
      my $dcz=$cz-$czz;
      my $dc=sqrt($dcx*$dcx+$dcy*$dcy+$dcz*$dcz);
      $dcx/=$dc;
      $dcy/=$dc;
      $dcz/=$dc;

      my $ddx=$cx-$cdx;
      my $ddy=$cy-$cdy;
      my $ddz=$cz-$cdz;
      my $dd=sqrt($ddx*$ddx+$ddy*$ddy+$ddz*$ddz);
      $ddx/=$dd;
      $ddy/=$dd;
      $ddz/=$dd;

      my $xx=$dcy*$ddz-$dcz*$ddy;
      my $xy=-$dcx*$ddz+$dcz*$ddx;
      my $xz=$dcx*$ddy-$dcy*$ddx;
      my $d=sqrt($xx*$xx+$xy*$xy+$xz*$xz);
      $xx/=$d;
      $xy/=$d;
      $xz/=$d;
 
      &checkIntersect($mol,$r,$cx,$cy,$cz,$xx,$xy,$xz);
    } elsif ($r->{name} eq "HIS" || $r->{name} eq "HSD" || $r->{name} eq "HSE") {
      for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
        my $n=$a->[$i]->{atomname};
        if ($n eq "CG" || $n eq "CE1" || $n eq "NE2" || $n eq "CD2" || $n eq "ND1") {
           $cx+=$a->[$i]->{xcoor};
           $cy+=$a->[$i]->{ycoor};
           $cz+=$a->[$i]->{zcoor};
           $nn++;
           if ($n eq "CG") {
             $czx=$a->[$i]->{xcoor};
             $czy=$a->[$i]->{ycoor};
             $czz=$a->[$i]->{zcoor};
           }
           if ($n eq "NE2") {
             $cdx=$a->[$i]->{xcoor};
             $cdy=$a->[$i]->{ycoor};
             $cdz=$a->[$i]->{zcoor};
           }
        } 
      }
      $cx/=$nn;
      $cy/=$nn;
      $cz/=$nn;

      my $dcx=$cx-$czx;
      my $dcy=$cy-$czy;
      my $dcz=$cz-$czz;
      my $dc=sqrt($dcx*$dcx+$dcy*$dcy+$dcz*$dcz);
      $dcx/=$dc;
      $dcy/=$dc;
      $dcz/=$dc;

      my $ddx=$cx-$cdx;
      my $ddy=$cy-$cdy;
      my $ddz=$cz-$cdz;
      my $dd=sqrt($ddx*$ddx+$ddy*$ddy+$ddz*$ddz);
      $ddx/=$dd;
      $ddy/=$dd;
      $ddz/=$dd;

      my $xx=$dcy*$ddz-$dcz*$ddy;
      my $xy=-$dcx*$ddz+$dcz*$ddx;
      my $xz=$dcx*$ddy-$dcy*$ddx;
      my $d=sqrt($xx*$xx+$xy*$xy+$xz*$xz);
      $xx/=$d;
      $xy/=$d;
      $xz/=$d;
 
      &checkIntersect($mol,$r,$cx,$cy,$cz,$xx,$xy,$xz);
    }  
  }
}

1;

sub checkIntersect  {
   my $mol=shift;
   my $r=shift;
   my $cx=shift;
   my $cy=shift;
   my $cz=shift;
   my $xx=shift;
   my $xy=shift;
   my $xz=shift;

#   printf STDERR "checking %d %s\n",$r->{num},$r->{name};
#   printf STDERR "  cofm %f %f %f\n",$cx,$cy,$cz;
#   printf STDERR "  normal %f %f %f\n",$xx,$xy,$xz;
    
   my $xd=sqrt($xx*$xx+$xy*$xy+$xz*$xz); 

   my @acand=();
   foreach my $tc ( @{$mol->{chain}} ) {
     foreach my $ta ( @{$tc->{atom}} ) {
       if ($ta->{resnum} ne $r->{num}) {
         my $dx=$cx-$ta->{xcoor};
         my $dy=$cy-$ta->{ycoor};
         my $dz=$cz-$ta->{zcoor};
         my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);
         if ($d<10) {
           push(@acand,$ta);
         }
       }
     }
   }
   for (my $ia=0; $ia<$#acand; $ia++) {
     for (my $ja=$ia+1; $ja<=$#acand; $ja++) {
       my $dx=$acand[$ia]->{xcoor}-$cx;
       my $dy=$acand[$ia]->{ycoor}-$cy;
       my $dz=$acand[$ia]->{zcoor}-$cz;
       my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);
       my $dotp1=($xx*$dx+$xy*$dy+$xz*$dz)/($xd*$d);       
       $dx=$acand[$ja]->{xcoor}-$cx;
       $dy=$acand[$ja]->{ycoor}-$cy;
       $dz=$acand[$ja]->{zcoor}-$cz;
       $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);
       my $dotp2=($xx*$dx+$xy*$dy+$xz*$dz)/($xd*$d);       
       
       $dx=$acand[$ia]->{xcoor}-$acand[$ja]->{xcoor};
       $dy=$acand[$ia]->{ycoor}-$acand[$ja]->{ycoor};
       $dz=$acand[$ia]->{zcoor}-$acand[$ja]->{zcoor};
       $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);

       if ($d<3.5 && abs($dotp1)>0.85 && abs($dotp2)>0.85 && $dotp1*$dotp2<0 && 
           abs($acand[$ia]->{resnum}-$acand[$ja]->{resnum})<=1 &&
           $acand[$ia]->{chain} eq $acand[$ja]->{chain}) {
           printf STDOUT " CROSS %s %d %s:: %s %d %s : %s %d %s : %f %f %f\n",$r->{seg},$r->{num},$r->{name},$acand[$ia]->{seg},$acand[$ia]->{resnum},$acand[$ia]->{atomname},$acand[$ja]->{seg},$acand[$ja]->{resnum},$acand[$ja]->{atomname},$dotp1,$dotp2,$d;
       }
     }
   }
}      

