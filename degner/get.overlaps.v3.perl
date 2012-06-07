#!/bin/perl

use strict;

my $snpfile = shift @ARGV;
my $genomepath = shift @ARGV;
my $readlength = shift @ARGV;
my $readquality = shift @ARGV;
my $seqfile;
my @chroms = (1..22, 'X', 'Y');
foreach my $chr (@chroms){
  read_files(($chr, $snpfile, $seqfile));
  print_overlaps($chr);
}

sub read_files {
  my $chr = shift @_;
  my $snpfile = shift @_;
  my $seqfile = shift @_;
  open(SNP, '<', $snpfile) or die "cannot open snpfile";
  open(SEQUENCE, "cat $genomepath/chr.$chr.fa |") or die "cannot open seqfile";
#DE: chr-formatted ref (UCSC):  open(SEQUENCE, "zcat $genomepath/chr$chr.fa.gz |") or die "cannot open seqfile";
}

sub print_overlaps {

  my $chr = shift @_;
  my @sequence_arr = <SEQUENCE>;
  chomp(@sequence_arr);
  shift(@sequence_arr); #remove fasta header
  my $sequence = join("", @sequence_arr);
#DE: UCSC formatted:  $sequence =~ s/>chr$chr//g;
  close(SEQUENCE);

  #DE: removes header: <SNP>;
  while (defined(my $line = <SNP>)) {
    chomp($line);
    my ($snpid, $thischr, $pos, $strand, $ref, $alt) = split(/ /, $line);

    my $plusbase;
#DE: chr-formatted ref:    if ($thischr eq "chr$chr") {
    if ($thischr eq $chr) {

      ###
      #Print short read identical to reference seq
      ###
      if ($strand eq "+") {
	$plusbase = $ref
      } elsif ($strand eq "-") {
	$plusbase = $ref;
	$plusbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
      }

      for (my $jbase=0; $jbase<$readlength; $jbase++) {
	my $short = substr($sequence, ($pos-($readlength-1-$jbase)-1), $readlength-1-$jbase) . $plusbase . substr($sequence, ($pos), $jbase);
	my $RC = $short;
	$RC =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
	$RC = reverse($RC);
	my $readname = join('_', ($snpid, $thischr, $pos, 'REF_+', $plusbase, $jbase));
	print('@', $readname, "\n", $short, "\n", '+', $readname, "\n", $readquality, "\n");
	$readname = join('_', ($snpid, $thischr, $pos, 'REF_-', $plusbase, $jbase));
	print('@', $readname, "\n", $RC, "\n", '+', $readname, "\n", $readquality, "\n");
#	print('@', "$snpid_$thischr_$pos_REF_+_$plusbase", "_$jbase", "\n$short\n", '+', "$snpid_$thischr_$pos_REF_+_$plusbase", "_$jbase", "\n$readquality\n");
#	print('@', "$snpid_$thischr_$pos_REF_-_$plusbase", "_$jbase", "\n$RC\n", '+', "$snpid_$thischr_$pos_REF_-_$plusbase", "_$jbase", "\n$readquality\n");
      }

      ###
      #Print short read with alternate allele
      ###
      if ($strand eq "+") {
	$plusbase = $alt
      } elsif ($strand eq "-") {
	$plusbase = $alt; 
	$plusbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
      }

      for (my $jbase=0; $jbase<$readlength; $jbase++) {
	my $short = substr($sequence, ($pos-($readlength-1-$jbase)-1), $readlength-1-$jbase) . $plusbase . substr($sequence, ($pos), $jbase);
	my $RC = $short;
	$RC =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
	$RC = reverse($RC);
	my $readname = join('_', ($snpid, $thischr, $pos, 'NONREF_+', $plusbase, $jbase));
	print('@', $readname, "\n", $short, "\n", '+', $readname, "\n", $readquality, "\n");
	$readname = join('_', ($snpid, $thischr, $pos, 'NONREF_-', $plusbase, $jbase)); 
	print('@', $readname, "\n", $RC, "\n", '+', $readname, "\n", $readquality, "\n");
#	print('@', "$snpid_$thischr_$pos_NONREF_+_$plusbase", "_$jbase", "\n$short\n", '+', "$snpid_$thischr_$pos_NONREF_+_$plusbase", "_$jbase", "\n$readquality\n");
#	print('@', "$snpid_$thischr_$pos_NONREF_-_$plusbase", "_$jbase", "\n$RC\n", '+', "$snpid_$thischr_$pos_NONREF_-_$plusbase", "_$jbase", "\n$readquality\n");
      }
    }
  }
  close(SNP);
}
