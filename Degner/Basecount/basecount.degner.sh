#Basecount (allelic count of each allele)
#DEP: dphetfilter.R
#DEP: mpileup.allelecount.pl
#DEP: vcfI162basecount.sh


########
#DEPTH and HETEROZYGOSITY filter
########

vcfdir='/proj/b2012046/rani/analysis/degner/varcall'

##
#VCF to tab format
##
PATH=${PATH}:/bubo/home/h26/edsgard/opt/vcftools_0.1.7/bin; export PATH
PERL5LIB=/bubo/home/h26/edsgard/opt/vcftools_0.1.7/lib; export PERL5LIB
vcftools --vcf allbams.vcf --extract-FORMAT-info DP --out genome
vcftools --vcf allbams.vcf --extract-FORMAT-info GT --out genome

##
#Filter
##
#OUT: *.het.pos
mindp=10
execdir='/bubo/home/h26/edsgard/glob/code/ase/sim'
Rscript ${execdir}/dphetfilter.R genome.DP.FORMAT genome.GT.FORMAT $mindp &
ls -1 *.het.pos | xargs -I% echo cat % '|' sed "'s/\./ /' >"%.hetvars >cmds.sh
sh cmds.sh
rename genome.DP.FORMAT.het.pos. . *.hetvars
rename .. . *.hetvars

#Dump "chr,pos,ref,alt" for these vars. (Used as input in snv.ase.R since the alt allele is seldom biallelic when running mpileup.allelecount.pl below.)
find $vcfdir -maxdepth 1 -name 'allbams.vcf'| xargs -I% echo "awk -F'\t' '{print \$1\".\"\$2, \$1, \$2, \$4, \$5;}' % | grep -v '^#' | sort -k 1,1 >"%.pos >cmds.sh
sh cmds.sh
ls -1 *.het.pos | xargs -I% echo sort % '>'%.sorted >cmds.sh
sh cmds.sh
ls -1 *.het.pos.sorted | sed 's/.genome.DP.FORMAT.het.pos.sorted//' >samples.list
cat samples.list | xargs -I% echo join -j1 %.genome.DP.FORMAT.het.pos.sorted allbams.vcf.pos '>' %.hetvars.alleles >cmds.sh
srun -p devel -t 1:00:00 -A b2012046 sh cmds.sh &
ls -1 *.hetvars.alleles | xargs -I% echo cut -d"' '" -f2-5 % "| tr ' ' '\t' >"%.tmp >cmds.sh
rename .hetvars.alleles.tmp .hetvars.alleles *.hetvars.alleles.tmp


#################################################
#BASECOUNT by MPILEUP (no varcall). 
#First four of I16 == DP4
#################################################
bamdir='/proj/b2012046/rani/analysis/degner/mergebam'
vcfdir='/proj/b2012046/rani/analysis/degner/varcall'
bcdir='/proj/b2012046/rani/analysis/degner/basecount'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/basecount'
scriptdir='/proj/b2012046/rani/scripts/degner/basecount'
projid='b2012046'
email='alva.rani@scilifelab.se'
time='10:00:00'
ref='/bubo/home/h24/alvaj/glob/annotation/degner/reference/Homo_sapiens.GRCh37.57.dna.concat.fa'

cd $scriptdir
find $bamdir -maxdepth 1 -name '*.bam' | xargs -I% basename % | sed 's/\.bam//' >samples.list

cd $scriptdir
find $bamdir -maxdepth 1 -name '*.bam' | xargs -I% basename % | sed 's/\.bam//' >samples.list
cat samples.list | xargs -I% echo 'perl' ${execdir}'/mpileup.allelecount.pl' ${bamdir}/%.bam ${vcfdir}/%.hetvars $ref $scriptdir $projid $email $time $bcdir >cmds.sh
sh cmds.sh
find $scriptdir -name '*mpileup.basecount*' | xargs -I% sbatch %

#Extract basecounts from vcf files (I16 field)

bcdir='/proj/b2012046/rani/analysis/degner/basecount'
execdir='/bubo/home/h26/edsgard/glob/code/ase'
find $bcdir -name '*nocall.vcf' >nocall.vcf.list
cat nocall.vcf.list | xargs -I% echo 'sh' ${execdir}/vcfI162basecount.sh % '>' %.basecount >cmds.sh
srun -p devel -t 1:00:00 -A b2012046 sh cmds.sh &

