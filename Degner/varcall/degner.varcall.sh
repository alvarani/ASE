#Varcalling with SAMtools
#DEP: mpileup_multisample.pl


###
#VarCall (SAMtools)
###

projid='b2012046'
email='alva.rani@scilifelab.se'
sbatchdir='/proj/b2012046/rani/scripts/degner/varscripts'
bamfile='/proj/b2012046/rani/analysis/degner/mergebam'
vcfdir='/proj/b2012046/rani/analysis/degner/varcall'
ref='/bubo/home/h24/alvaj/glob/annotation/degner/reference/Homo_sapiens.GRCh37.57.dna.concat.fa'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/varcall'

cd $sbatchdir
# need to redo it again with merged bam files
#create file listing the bam files
allbamsfile=${vcfdir}/allbams.list
find $bamfile -maxdepth 1 -name '*.bam' >$allbamsfile

#create region files
module load bioinfo-tools
module load samtools/0.1.18
srun -p devel -A b2012046 -t 1:00:00 samtools faidx $ref &
vcfutils.pl splitchr -l 6600000 ${ref}.fai >${vcfdir}/genome.480.regs


#Run mpileup
cat ${vcfdir}/genome.480.regs | xargs -I% echo 'perl' ${execdir}/mpileup_multisample.pl $allbamsfile % $sbatchdir $projid $email $ref >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.mpileup.sbatch' | xargs -I% sbatch %


#check status:
cd $vcfdir
lt *.stderr >allerrfiles.tmp
squeue -u alvaj | grep -v JOBID | awk '{print $1;}' | xargs -I% grep % ${vcfdir}/info/allerrfiles.tmp

#Catenate all region-specific vcfs
find $vcfdir -name 'allbams.*.vcf' | xargs cat | grep -v '^#' >${vcfdir}/allbams.vcf.tmp
cat 'allbams.list.1:1-6600000.vcf' | grep '^#' >header.tmp
cat header.tmp allbams.vcf.tmp >allbams.vcf
rm header.tmp allbams.vcf.tmp
wc -l allbams.vcf 
#298942 

#Rm slashes from vcf sample header
cd ${vcfdir}
cat allbams.vcf | grep '^#' >header.tmp #/proj/b2012046/rani/analysis/degner/mergebam
cat header.tmp | sed 's/\/proj\/b2012046\/rani\/analysis\/degner\/mergebam\///g' | sed 's/\.bam//g' >file.tmp
mv file.tmp header.tmp
grep -v '^#' allbams.*.vcf | sed 's/^.*\.vcf://' >allbams.vcf.noheader
cat header.tmp allbams.vcf.noheader >allbams.vcf
rm allbams.vcf.noheader header.tmp
