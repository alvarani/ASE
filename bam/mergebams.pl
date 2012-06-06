#!/usr/bin/perl

#find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -name 'accepted_hits.bam' | awk -F'/' '{print "perl /bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/rmdups.pl", $0, $10, "/proj/b2011075/analysis/hs_pipe/outscripts/run1 b2010035 daniel.edsgard@scilifelab.se";}' >rmdups.sh

use strict;

#my ($sampleid, $bamdir, $outdir, $sbatch_dir, $aid, $em) = @ARGV;

mergesams(@ARGV);

sub mergesams {

  use File::Spec::Functions;
  use File::Basename;

  my $sampleid = shift @_;
  my $bamdir = shift @_;
  my $outdata_dir = shift @_;
  my $sbatch_dir = shift @_;
  my $aid = shift @_;
  my $em = shift @_;
  my $time = shift @_;

  if(!defined($time)){
    $time = '23:59:00';
  }

  my $sbatch_file = catfile($sbatch_dir, "$sampleid.mergesams.sbatch");
  my $outdata_infodir = catfile($outdata_dir, 'info');

  #get bamfiles to be merged
  my @bamfiles = `find $bamdir -name '*.bam' | grep $sampleid`;
  chomp(@bamfiles);

  #set outfile
  my $bamfile_merged = catfile($outdata_dir, $sampleid . '.bam');

  #make dirs
  `mkdir -p $outdata_infodir;`; #Creates the data and info dir
  `mkdir -p $sbatch_dir;`; #Creates the sbatch-script dir

  open (SBATCH, '>', $sbatch_file) or die("Can't write to $sbatch_file: $!\n");

  #print sbatch vars
  print SBATCH "#!/bin/bash -l", "\n";
  print SBATCH "#SBATCH -A ", $aid, "\n";
  print SBATCH "#SBATCH -t ", $time, "\n";
  print SBATCH "#SBATCH -J mergesams_", $sampleid, "\n";
  print SBATCH "#SBATCH -p core -n 1", "\n";
  print SBATCH "#SBATCH -e $outdata_infodir/$sampleid" . '.mergesams.jid_%.stderr' . "\n";
  print SBATCH "#SBATCH -o $outdata_infodir/$sampleid" . '.mergesams.jid_%.stdout' . "\n";
  unless ($em eq 0) {	
    print SBATCH "#SBATCH --mail-type=All", "\n";
    print SBATCH "#SBATCH --mail-user=$em", "\n\n";	
  }

  #print vars    
  print SBATCH 'module load bioinfo-tools' . "\n";
  print SBATCH 'module load picard/1.41' . "\n";

  #print cmd
  print SBATCH 'java -Xmx2g -jar /bubo/sw/apps/bioinfo/picard/1.41/MergeSamFiles.jar INPUT=' . join(' INPUT=', @bamfiles) . ' OUTPUT=' . $bamfile_merged . ' USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT' . "\n";
  close(SBATCH);        
}
