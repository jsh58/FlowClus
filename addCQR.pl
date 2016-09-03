#!/usr/bin/perl

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Dec. 2013

# FlowClus can be used to denoise Ion Torrent data,
#   but it requires an sff.txt file as an input.
# The sff.txt file produced by Ion Torrent is
#   similar to that of 454, but not identical. In
#   particular, the "Clip Qual Right" data field
#   is not set. Unfortunately, FlowClus uses this
#   number when retrieving the flowgram.
# This script creates a new sff.txt file with
#   updated "Clip Qual Right" values, based on the
#   "# of Bases" field. FlowClus can then be used
#   to denoise the new sff.txt file.

use strict;
use warnings;

die "Usage:  perl addCQR.pl <inputSFF> <outputSFF>\n",
    "Error!  Need input and output sff.txt files\n"
  if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

open(IN, $ARGV[0]) || die "Error!  Cannot open $ARGV[0]\n$!\n";
die "Error!  $ARGV[1] already exists\n" if (-f $ARGV[1]);
open(OUT, ">$ARGV[1]");

my $x = 0;
while (my $line = <IN>) {
  chomp $line;
  my @spl = split(":", $line);
  if (scalar @spl > 1 && $spl[0] eq "  Clip Qual Right") {
    die "No value given\n" if (! $x); 
    print OUT "$spl[0]:$x\n";
    $x = 0;
  } else {
    if (scalar @spl > 1 && $spl[0] eq "  # of Bases") {
      $x = substr($spl[1], 5);
    }
    print OUT "$line\n";
  }
}

close IN;
close OUT;
