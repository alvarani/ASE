#!/usr/bin/perl


#-q1 mapQ filter to get "unique" alignments.

use strict;
use File::Basename;
use File::Spec;

my ($bamfile, $varsfile, $ref, $ods, $aid, $em, $time, $outdir) = @ARGV;

if(!(defined($time))){
  $time = '2:00:00';
}
if(!(defined($outdir))){
  $outdir = dirname($bamfile);
}


mpileup(($bamfile, $varsfile, $ref, $ods, $aid, $em, $time, $outdir));

sub mpileup {

  use File::Spec::Functions;
  use File::Basename;

  my $bamfile = shift @_;
  my $varsfile = shift @_;
  my $ref = shift @_;
  my $ods = shift @_;
  my $aid = shift @_;
  my $em = shift @_;
  my $time = shift @_;
  my $outdir = shift @_;

  my $bamfile_base = basename($bamfile);
  my $sbatch_base = $bamfile;
  $sbatch_base =~ s/\///;
  $sbatch_base =~ s/\//\./g;
  my $sbatch_dir = $ods;
  my $sbatch_file = catfile($sbatch_dir, $sbatch_base . '.mpileup.basecount.sbatch');
  my $outdata_infodir = catfile($outdir, 'info');

#  print $sbatch_file . "\n";

  #make dirs
  `mkdir -p $outdata_infodir;`; #Creates the data and info dir
  `mkdir -p $sbatch_dir;`; #Creates the sbatch-script dir

  open (SBATCH, '>', $sbatch_file) or die("Can't write to $sbatch_file: $!\n");

  #print sbatch varsfile
  print SBATCH "#!/bin/bash -l", "\n";
  print SBATCH "#SBATCH -A ", $aid, "\n";
  print SBATCH "#SBATCH -t ", $time, "\n";
  print SBATCH "#SBATCH -J multimpileup", "\n";
  print SBATCH "#SBATCH -p core -n 1", "\n";
  print SBATCH "#SBATCH -e $outdata_infodir/$sbatch_base.mpileup.basecount" . '.jid_%j.stderr', "\n";
  print SBATCH "#SBATCH -o $outdata_infodir/$sbatch_base.mpileup.basecount" . '.jid_%j.stdout', "\n";
  unless ($em eq 0) {	
    print SBATCH "#SBATCH --mail-type=All", "\n";
    print SBATCH "#SBATCH --mail-user=$em", "\n\n";	
  }

  #print varsfile    
  print SBATCH 'module load bioinfo-tools' . "\n";
  print SBATCH 'module load samtools/0.1.18' . "\n";

  #print cmd
  my $basecountfile = catfile($outdir, $bamfile_base . '.mpileup.nocall.vcf');

  print SBATCH 'samtools mpileup -q1 -d10000 -L10000 -DSugf ' . $ref . ' -l ' . $varsfile . ' ' . $bamfile . q( | bcftools view - >) . $basecountfile . "\n";  

  close(SBATCH);        
}
