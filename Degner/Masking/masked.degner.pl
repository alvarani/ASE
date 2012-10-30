
# perl masked.degner.pl /proj/b2010035/edsgard/annotations/1000g/2010_nov/EUR.snps.2010_11.tab.nohomo.trifixed /bubo/nobackup/uppnex/reference/Homo_sapiens/GRCh37/concat/Homo_sapiens.GRCh37.57.dna.concat.fa /bubo/home/h24/alvaj/glob/annotation/degner/reference/Homo_sapiens.maskedRefe.GRCh37.57.fa


#!/bin/perl
use strict;
use warnings;
no warnings ('uninitialized', 'substr');


my $snpfile = $ARGV[0];
my $sequencefile = $ARGV[1];
my $chromosomename=$ARGV[2];
#my $chr;

open(SNP, "$snpfile");
open(SEQUENCE,"genome.fa");
open(NRSEQ, ">$chromosomename");

print_chromosome();

sub print_chromosome {
my @line;
my $line;
my $plusbase;
my $otherbase;
my $test_string;
my $NRsequence;
my $chr; my $ky;
my $seq; my $chr_snp;
my %seq;
#<SEQUENCE>;
while (<SEQUENCE>){
chomp;
$NRsequence = $_;
 $chr = $1 if ($NRsequence =~ />(\d+)\s(.*)/) || ($NRsequence =~ />(\w+)\s(.*)/ );
$seq{$chr}.=$NRsequence unless ($NRsequence =~ />/);  
}

#close(SEQUENCE);
my $a = keys %seq;
print "keys:",$a;
#print $_,"\n\n" for keys %seq;

#<SNP>;
while($line=<SNP>){
@line = split(/\t/, $line);
#print @line, "\n";
#if($line[3] eq "+"){
$chr_snp =$line[0];
$plusbase=$line[4]; # reference
$otherbase=$line[3]; # alternative

print "Ref:$plusbase,Alt:$otherbase\n";

#}
#elsif($line[3] eq "-"){
#$plusbase = $line[4];
#$otherbase=$line[5];

#$plusbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
#$otherbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
#}

#print "$line[2]\n"; #
print(substr($seq{$chr_snp}, $line[2]-5, 10), "\n");
$test_string = substr($seq{$chr_snp}, $line[2]-1, 1);
if( uc($plusbase) eq uc($test_string) ) {
print "OK\n";
if(uc($plusbase) ne "C" && uc($otherbase) ne "C"){
substr $seq{$chr_snp}, $line[2]-1, 1, "C";
} elsif(uc($plusbase) ne "G" && uc($otherbase) ne "G"){
substr $seq{$chr_snp}, $line[2]-1, 1, "G";
} elsif(uc($plusbase) ne "A" && uc($otherbase) ne "A"){
substr $seq{$chr_snp}, $line[2]-1, 1, "A";
} elsif(uc($plusbase) ne "T" && uc($otherbase) ne "T"){
substr $seq{$chr_snp}, $line[2]-1, 1, "T";
}
} elsif($test_string eq "N"){
print "Base was N in reference\n";
}
else {
die "THE BASE GIVEN WAS NOT THE REFERENCE BASE, ABORTING MISSION!";
print "THE BASE GIVEN WAS NOT THE REFERENCE BASE, ABORTING MISSION!";
}
print(substr($seq{$chr_snp}, $line[2]-5, 10),"\n\n");
}
#print NRSEQ (">chr\n$seq{$chr_snp}\n");
#print NRSEQ ">chr",$_,"\n",$seq{$_},"\n" for each %seq;
foreach $ky (sort keys %seq) {print NRSEQ ">chr$ky\n$seq{$ky}\n";}
}

