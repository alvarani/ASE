#!/usr/bin/perl



#EX: perl rm.multimappedreads.pl *.sam *.nonuniq



use strict;



my $samfile = shift @ARGV;

my $rmfile = shift @ARGV;



#slurp rmfile

my %reads2rm;

open(RMFILE, '<', $rmfile) or die("Couldn't open $rmfile, $!");

while(defined(my $line = <RMFILE>)){

  chomp($line);

  my $read = $line;

  $reads2rm{$read} = ();

}

close(RMFILE);



open(SAM, '<', $samfile);

while(defined(my $line = <SAM>)){

  chomp($line);

  my @cols = split(/\t/, $line);

  my $read = $cols[0];


