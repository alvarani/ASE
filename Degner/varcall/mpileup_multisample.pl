#!/usr/bin/perl

#find '/proj/b2011075/analysis/hs_pipe/outdata/run1/varcalls' -name 'allbams.chr_*' | xargs -I% echo 'perl /bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/mpileup_multisample.pl' % '/proj/b2011075/analysis/hs_pipe/outscripts/run1/varcalls b2010035 daniel.edsgard@scilifelab.se' >cmds.sh

#-q1 mapQ filter to get "unique" alignments.
#No max dp varFilter since RNA-seq data.

use strict;

my ($bamfilelist, $region, $ods, $aid, $em, $ref, $time, $partition) = @ARGV;

if(!(defined($ref))){
  $ref = '/bubo/home/h26/edsgard/glob/annotation/human/Homo_sapiens.GRCh37.57.dna.concat.fa';
}
if(!(defined($time))){
  $time = '23:00:00';
}
if(!(defined($partition))){
  $partition = 'node';
}

mpileup(($bamfilelist, $region, $ods, $aid, $em, $ref, $time, $partition));

sub mpileup {

  use File::Spec::Functions;
  use File::Basename;

  my $bamfilelist = shift @_;
  my $region = shift @_;
  my $ods = shift @_;
  my $aid = shift @_;
  my $em = shift @_;
  my $ref = shift @_;
  my $time = shift @_;
  my $partition = shift @_;

  my $bamfiles_base = basename($bamfilelist);
  my $sbatch_dir = $ods;
  my $sbatch_file = catfile($sbatch_dir, $bamfiles_base . '.' . $region . '.mpileup.sbatch');
  my $outdata_dir = dirname($bamfilelist);
  my $outdata_infodir = catfile($outdata_dir, 'info');

  #make dirs
  `mkdir -p $outdata_infodir;`; #Creates the data and info dir
  `mkdir -p $sbatch_dir;`; #Creates the sbatch-script dir

  open (SBATCH, '>', $sbatch_file) or die("Can't write to $sbatch_file: $!\n");

  #print sbatch vars
  print SBATCH "#!/bin/bash -l", "\n";
  print SBATCH "#SBATCH -A ", $aid, "\n";
  print SBATCH "#SBATCH -t ", $time, "\n";
  print SBATCH "#SBATCH -J multimpileup", "\n";
  print SBATCH '#SBATCH -p ', $partition, ' -n 1', "\n";
  print SBATCH "#SBATCH -e $outdata_infodir/$bamfiles_base.$region.multimpileup" . '.jid_%j.stderr', "\n";
  print SBATCH "#SBATCH -o $outdata_infodir/$bamfiles_base.$region.multimpileup" . '.jid_%j.stdout', "\n";
  unless ($em eq 0) {	
    print SBATCH "#SBATCH --mail-type=All", "\n";
    print SBATCH "#SBATCH --mail-user=$em", "\n\n";	
  }

  #print vars    
  print SBATCH 'module load bioinfo-tools' . "\n";
  print SBATCH 'module load samtools/0.1.18' . "\n";
  print SBATCH 'echo "Running on: $(hostname)"',"\n\n";

  #print cmd
  my $vcffile = $bamfilelist . '.' . $region . '.vcf';
  print SBATCH 'samtools mpileup -q1 -d10000 -L10000 -DSugf ' . $ref . ' -r ' . $region . ' -b ' . $bamfilelist. ' | bcftools view -vcg - | vcfutils.pl varFilter >' . $vcffile . "\n";

  close(SBATCH);        
}
