#Varcalling with SAMtools
#DEP: mpileup_multisample.pl


###
#VarCall (SAMtools)
###

projid='b2012046'
email='alva.rani@scilifelab.se'
sbatchdir='/proj/b2012046/rani/scripts/varcall'
bamdir='/proj/b2012046/rani/analysis/gsnap'
vcfdir='/proj/b2012046/edsgard/ase/sim/data/varcalls'
ref='/bubo/home/h26/edsgard/glob/annotation/human/Homo_sapiens.GRCh37.57.dna.concat.fa'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/varcall'

cd $sbatchdir

#create file listing the bam files
allbamsfile=${vcfdir}/allbams.list
find $bamdir -maxdepth 1 -name '*.bam' >$allbamsfile

#create region files
module load bioinfo-tools
module load samtools/0.1.18
srun -p devel -A b2010035 -t 1:00:00 samtools faidx $ref &
vcfutils.pl splitchr -l 6600000 ${ref}.fai >${vcfdir}/genome.480.regs

#Run mpileup
cat ${vcfdir}/genome.480.regs | xargs -I% echo 'perl' ${execdir}/mpileup_multisample.pl $allbamsfile % $sbatchdir $projid $email $ref >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.mpileup.sbatch' | xargs -I% sbatch %

#Check status:
cd $vcfdir
lt *.stderr >allerrfiles.tmp
squeue -u edsgard | grep -v JOBID | awk '{print $1;}' | xargs -I% grep % ${vcfdir}/info/allerrfiles.tmp

#Catenate all region-specific vcfs
find $vcfdir -name 'allbams.*.vcf' | xargs cat | grep -v '^#' >${vcfdir}/allbams.vcf.tmp
cat allbams.list.MT:1-16569.vcf | grep '^#' >header.tmp
cat header.tmp allbams.vcf.tmp >allbams.vcf
rm header.tmp allbams.vcf.tmp
wc -l allbams.vcf #145615

#Rm slashes in vcf sample header
vcfdir='/proj/b2012046/edsgard/ase/sim/data/varcalls'
cd ${vcfdir}
cat allbams.vcf | grep '^#' >header.tmp
cat header.tmp | sed 's/\/proj\/b2012046\/edsgard\/ase\/sim\/data\/tophat\///g' | sed 's/\.bam//g' >file.tmp
mv file.tmp header.tmp
grep -v '^#' allbams.vcf >allbams.vcf.noheader
cat header.tmp allbams.vcf.noheader >allbams.vcf
rm allbams.vcf.noheader header.tmp
