#!/usr/bin/perl

use strict;
use File::Basename;
use File::Spec::Functions;


my ($sample, $fastqprefix, $fastqsuffix, $outdir, $threads, $gtf, $ref, $isizefile, $ods, $projid, $email, $time, $readlen, $isizedev) = @ARGV;

#set vars
my $fastqfile_1 = $fastqprefix . '1' . $fastqsuffix;
#my $fastqfile_2 = $fastqprefix . '2' . $fastqsuffix;
my $fastqprefix_base = basename($fastqprefix);
my $outdir = catfile($outdir, $fastqprefix_base);
$outdir =~ s/\.$//;

#set sbatch file
my $sbatch_file = catfile($ods, $fastqprefix_base . $fastqsuffix . '.tophat.sbatch');
$sbatch_file =~ s/\.\./\./;

#set insert size (NB: actually inner distance)
my $isize;
open(ISIZE, '<', $isizefile) or die("Couldnt open $isizefile, $!\n");
while(defined(my $line = <ISIZE>)){
  chomp($line);
#  print 'DEBUG: ' . $line . "\n";
  my ($jsample, $isize_jsample) = split(/\t/, $line);
#  print 'DEBUG: jsample: ' . $jsample . "\n";
#  print 'DEBUG: sample: ' . $sample . "\n";
#  my ($isize_jsample, $isize_jsample, $isizedev_jsample) = split(/ /, $line);
  if($jsample =~ $sample){
    $isize = $isize_jsample;
    $isize = $isize - 2 * $readlen; #Tophat uses the inner distance: isize - 2*readlen
#    $isizedev = $isizedev_jsample;
  }
}
close(ISIZE);

#create outdirs
my $cmd = 'mkdir -p ' . catfile($outdir, 'log');
system($cmd);
my $cmd = 'mkdir -p ' . $ods;
system($cmd);

#print sbatch file
open(SBATCH, '>', $sbatch_file) or die("Couldn't open $sbatch_file");
print SBATCH '#!/bin/bash -l' . "\n";
print SBATCH '#SBATCH -A ' . $projid . "\n";
print SBATCH '#SBATCH -t ' . $time . "\n";
print SBATCH '#SBATCH -J tophat' . "\n";
print SBATCH '#SBATCH -p node -n ' . $threads . "\n";
print SBATCH '#SBATCH -e ' . $outdir . '/log/tophat.jid_%j.stderr' . "\n";
print SBATCH '#SBATCH -o ' . $outdir . '/log/tophat.jid_%j.stdout' . "\n";
print SBATCH '#SBATCH --mail-type=All' . "\n";
print SBATCH '#SBATCH --mail-user=' . $email . "\n";

#set vars
print SBATCH 'module load bioinfo-tools' . "\n";
print SBATCH 'module load samtools/0.1.18' . "\n";
print SBATCH 'module load tophat/1.4.0' . "\n";

#print cmd
print SBATCH 'tophat --solexa1.3-quals -p ' . $threads . ' -o ' . $outdir . ' --GTF ' . $gtf . " -r '" . $isize . "' --mate-std-dev " . $isizedev . ' ' . $ref . ' ' . $fastqfile_1 . ' ' . "\n";
close(SBATCH);

