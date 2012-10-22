#!/usr/bin/perl

#perl rmdups.pl /proj/b2011075/analysis/hs_pipe/outdata/run1/1_unstim_12h/tophat/1_unstim_12h.lane_1.index_8.readp_1.fastq.1.filter.fastq/accepted_hits.bam 1_unstim_12h.lane_1.index_8 '/proj/b2011075/analysis/hs_pipe/outscripts/run1' 'b2010035' 'daniel.edsgard@scilifelab.se'
#find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -name 'accepted_hits.bam' | awk -F'/' '{print "perl /bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/rmdups.pl", $0, $10, "/proj/b2011075/analysis/hs_pipe/outscripts/run1 b2010035 daniel.edsgard@scilifelab.se";}' >rmdups.sh

use strict;

#my ($bamfile, $sampleid, $ods, $aid, $em) = @ARGV;

rmDups(@ARGV);

sub rmDups {

  use File::Spec::Functions;
  use File::Basename;

  my $bamfile = shift @_;
  my $sampleid = shift @_;
  my $ods = shift @_;
  my $aid = shift @_;
  my $em = shift @_;
  my $time = '9:00:00';

  my $sbatch_dir = catfile($ods, $sampleid);
  my $sbatch_file = catfile($sbatch_dir, "$sampleid.rmdups.sbatch");
  my $outdata_dir = dirname($bamfile);
  my $outdata_infodir = catfile($outdata_dir, 'log');

  #make dirs
  `mkdir -p $outdata_infodir;`; #Creates the data and info dir
  `mkdir -p $sbatch_dir;`; #Creates the sbatch-script dir

  open (SBATCH, '>', $sbatch_file) or die("Can't write to $sbatch_file: $!\n");

  #print sbatch vars
  print SBATCH "#!/bin/bash -l", "\n";
  print SBATCH "#SBATCH -A ", $aid, "\n";
  print SBATCH "#SBATCH -t ", $time, "\n";
  print SBATCH "#SBATCH -J rmDups_", $sampleid, "\n";
  print SBATCH "#SBATCH -p core -n 1", "\n";
  print SBATCH "#SBATCH -e $outdata_infodir/$sampleid.rmdups.stderr", "\n";
  print SBATCH "#SBATCH -o $outdata_infodir/$sampleid.rmdups.stdout", "\n";
  unless ($em eq 0) {	
    print SBATCH "#SBATCH --mail-type=All", "\n";
    print SBATCH "#SBATCH --mail-user=$em", "\n\n";	
  }

  #print vars    
  print SBATCH 'module load bioinfo-tools' . "\n";
  print SBATCH 'module load picard/1.41' . "\n";
  print SBATCH 'echo "Running on: $(hostname)"',"\n\n";

  #print cmd
  my $bamfile_sorted = $bamfile . '.sorted';
  my $bamfile_nodup = $bamfile_sorted . '.nodup';
  print SBATCH 'java -Xmx2g -jar /bubo/sw/apps/bioinfo/picard/1.41/SortSam.jar INPUT=' . $bamfile . ' OUTPUT=' . $bamfile_sorted . ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT' . "\n";
  print SBATCH 'java -Xmx2g -jar /bubo/sw/apps/bioinfo/picard/1.41/MarkDuplicates.jar INPUT=' . $bamfile_sorted . ' OUTPUT=' . $bamfile_nodup . ' ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=' . $sampleid . '_picardDup_metrics VALIDATION_STRINGENCY=LENIENT' . "\n";

  close(SBATCH);        
}
