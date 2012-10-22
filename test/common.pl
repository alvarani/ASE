#!/usr/bin/perl
$\=undef;
use strict;
use warnings;
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
open (FILE1, $file1) || die "Could not open file \n";
open (FILE2, $file2) || die "Could not open file \n";
my $outfile = $ARGV[2];
my @outlines;
my @outlines2;
my @outlines_cor;
my @outlines_swiss;
my %hash = ();
my %hash_2 = ();
open (OUTFILE, ">$outfile") or die "Cannot open $outfile for writing \n";
#x.txt
while (<FILE1>) {
    my @array = split(/\t/);
    $hash{$array[0]} = join("\t",@array);
}
 
#y.txt
while (<FILE2>) {
    #next unless (/^[A-Z]/);#skip lines that do not start with an uppercase alpha
    chomp;
    my $col1 = (split(/\t/))[0]; 
    if (exists ($hash{$col1})) {
        print OUTFILE $hash{$col1};
    }else {#print OUTFILE $col1."\n";


}
}
close OUTFILE;
close FILE1;
close FILE2;
