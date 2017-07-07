#!/usr/bin/env perl

while (<>) {
  if (/ATOM/) {
    chomp;
    my $atomname=substr($_,12,4);
    $atomname=~s/ //g;
    $atomname=~s/^\d//g;
    if ($atomname eq "CLA") {
      $element="CL";
    } elsif ($atomname eq "SOD") {
      $element="NA";
    } elsif ($atomname eq "CAL") {
      $element="CA";
    } else {
      $element=" ".substr($atomname,0,1);
    }
    my $line=" "x78;
    substr($line,0,length($_))=$_;
    substr($line,76,2)=$element; 
    print $line,"\n";
  } else {
    print;
  }
}
      
