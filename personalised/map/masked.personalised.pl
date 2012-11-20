
# perl masked.personalised.pl /proj/b2011075/analysis/sim/data/synt/snps/synt.snps.tab.1 /bubo/nobackup/uppnex/reference/Homo_sapiens/GRCh37/concat/Homo_sapiens.GRCh37.57.dna.concat.fa /bubo/home/h24/alvaj/glob/annotation/degner/personalised/Homo_sapiens.maskedRefe.personalised.GRCh37.57.fa



#!/bin/perl
use strict;
use warnings;
no warnings ('uninitialized', 'substr');


my $snpfile = $ARGV[0];
my $sequencefile = $ARGV[1];
my $chromosomename=$ARGV[2];
#my $chr;

open(SNP, "$snpfile");
open(SEQUENCE,"$sequencefile");
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
my $chrnum; my $local_seq;
#<SEQUENCE>;
while (<SEQUENCE>){
chomp;
$NRsequence = $_;
#print $_,"\n";
 $chr = $1 if ($NRsequence =~ />(\d+)\s(.*)/) || ($NRsequence =~ />(\w+)\s(.*)/ );
#print "$chr\n";
$seq{$chr}.=$NRsequence unless ($NRsequence =~ />/);
}

close(SEQUENCE);
print "$seq{1}\n";
my $a = keys %seq;
print "$a keys: \n";


#<SNP>;
while($line=<SNP>){
chomp;
@line = split(/\t/, $line);

#if($line[3] eq "+"){
$chr_snp =$line[2]; # snp position
$chrnum = $line[1];
$plusbase=$line[3]; #reference
$otherbase=$line[4]; # alternative
my $local_seq = $seq{$chrnum};

#print "$line[2]\n"; #
print(substr($local_seq, $line[2]-5, 10), "\n"); 
$test_string = substr($local_seq, $line[2]-1, 1);
if(uc($plusbase) eq uc($test_string)){
print "OK\n";
if(uc($plusbase) ne "C" && uc($otherbase) ne "C"){
substr $local_seq, $line[2]-1, 1, "C";
} elsif(uc($plusbase) ne "G" && uc($otherbase) ne "G"){
substr $local_seq, $line[2]-1, 1, "G";
} elsif(uc($plusbase) ne "A" && uc($otherbase) ne "A"){
substr $local_seq, $line[2]-1, 1, "A";
} elsif(uc($plusbase) ne "T" && uc($otherbase) ne "T"){
substr $local_seq, $line[2]-1, 1, "T";
}
} elsif($test_string eq "N"){
print "Base was N in reference\n";
}
else {
die "THE BASE GIVEN WAS NOT THE REFERENCE BASE, ABORTING MISSION!";
#print "THE BASE GIVEN WAS NOT THE REFERENCE BASE, ABORTING MISSION!";
}
print(substr($local_seq, $line[2]-5, 10),"\n\n");
}
#print NRSEQ (">chr\n$seq{$chr_snp}\n");
#print NRSEQ ">chr",$_,"\n",$seq{$_},"\n" for each %seq;
foreach $ky (sort keys %seq) {print NRSEQ ">chr$ky\n$seq{$ky}\n";}
}


