#!/usr/bin/perl

#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" | sed 's/^chr//' | sed 's/^M/MT/' > $chrsizes
#prg='/bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/coverage.pl'
#find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -name 'merged.bam' | awk -F'/' -v ods="$ods" -v projid="$projid" -v em="$em" -v prg="$prg" '{print "perl", prg, $0, $8, ods, projid, em;}' >bedtoolscov.sh

use strict;

#my ($bamfile, $sampleid, $ods, $aid, $em) = @ARGV;

my $bed='/bubo/home/h26/edsgard/glob/annotation/human/CCDS-hg19.chrstripped.bed';
my $gtf='/bubo/home/h26/edsgard/glob/annotation/human/Homo_sapiens.GRCh37.59.gtf';
my $genome='/bubo/home/h26/edsgard/glob/annotation/human/hg19.chrsizes.chrstripped';

coverage((@ARGV, $bed, $gtf, $genome));

sub coverage {

  use File::Spec::Functions;
  use File::Basename;

  my $bamfile = shift @_;
  my $sampleid = shift @_;
  my $ods = shift @_;
  my $em = shift @_;
  my $bed = shift @_;
  my $gtf = shift @_;
  my $genome = shift @_;
  my $time = '3:00:00';

  my $sbatch_dir = catfile($ods, $sampleid, 'coverage');
  my $sbatch_file = catfile($sbatch_dir, "$sampleid.coverage.sbatch2");
  my $outdata_dir = catfile(dirname(dirname($bamfile)), 'coverage');
  my $outdata_infodir = catfile($outdata_dir, 'info');

  #make dirs
  `mkdir -p $outdata_infodir;`; #Creates the data and info dir
  `mkdir -p $sbatch_dir;`;

  open (SBATCH, '>', $sbatch_file) or die("Can't write to $sbatch_file: $!\n");

  #print sbatch vars
 print SBATCH "#!/bin/bash -l", "\n";
  print SBATCH "#SBATCH -A ", $ods, "\n";
  print SBATCH "#SBATCH -t ", $time, "\n";
  print SBATCH "#SBATCH -J Cov_", $sampleid, "\n";
  print SBATCH "#SBATCH -p node -n 1", "\n";
  print SBATCH "#SBATCH -e $outdata_infodir/$sampleid.bedtools.cov.stderr2", "\n";
  print SBATCH "#SBATCH -o $outdata_infodir/$sampleid.bedtools.cov.stdout2", "\n";

  unless ($em eq 0) {	
    print SBATCH "#SBATCH --mail-type=All", "\n";
    print SBATCH "#SBATCH --mail-user=$em", "\n\n";	
  }

  #print vars    
  print SBATCH 'module load bioinfo-tools' . "\n";
  print SBATCH 'module load BEDTools/2.11.2' . "\n";
  print SBATCH 'echo "Running on: $(hostname)"',"\n\n";

  #print cmd
  print SBATCH 'coverageBed -abam ' . $bamfile . ' -b ' . $bed . ' -hist -split >' . catfile($outdata_dir, 'ccds.bedtools.out2') . "\n";
  print SBATCH 'coverageBed -abam ' . $bamfile . ' -b ' . $gtf . ' -hist -split >' . catfile($outdata_dir, 'ensembl.bedtools.out2') . "\n";
  print SBATCH 'genomeCoverageBed -ibam ' . $bamfile . ' -g ' . $genome . ' -split -max 1000 >' . catfile($outdata_dir, 'genome.bedtools.out2') . "\n";

  close(SBATCH);        
}
