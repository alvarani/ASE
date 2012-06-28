#! /bin/bash -l
#Map with tophat
#DEP: tophat.gensbatch.pl


###
#Tophat
###
#Longest run: 5h38min, 6G file for each readpair: 8_LPS.readp_1.2.ill.fastq
email='alva.rani@scilifelab.se'
projid='b2012046'
#execdir='/bubo/home/h24/alvaj/glob/code/ASE/map'
fastqdir='/proj/b2012046/edsgard/ase/sim/data/degner'
#fastqdir='/proj/b2012046/rani/data/fastq'
#fastqfile='testsingle.bqB.fastq'
outdir='/proj/b2012046/rani/analysis/monte'
sbatchdir='/proj/b2012046/rani/scripts/monte'

#Vars.
#threads='8'
annotdir='/bubo/home/h26/edsgard/glob/annotation/human/'
gtf=${annotdir}/'Homo_sapiens.GRCh37.59.gtf'
ref='/bubo/home/h26/edsgard/glob/annotation/bowtie/concat.fa'
#time='01:00:00'

#Longest run: 7h10min: Input: ~7.3G per readpair fastq file

#readlen='100'
#isize='175'
#isizedev='75' 

#very rough approx from the bioanalyzer plots. Also, we use trimming which further adds to the deviation.

cd $sbatchdir
#rm cmds.sh
#samples=(`cat samples.list`)

#for sample in ${samples[@]}
#do
#find $fastqdir -maxdepth 1 -name ${samples}'*' | sort -u >fastqfiles.list

(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2012046
#SBATCH -t 01:00:00
#SBATCH -J tophat
#SBATCH -p node -n 8
#SBATCH -e /proj/b2012046/rani/analysis/monte/log/tophat.jid_%j.stderr
#SBATCH -o /proj/b2012046/rani/analysis/monte/log/tophat.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=alva.rani@scilifelab.se
module load bioinfo-tools
module load samtools/0.1.18
module load tophat/1.4.0
tophat --solexa1.3-quals -G $gtf -o  $outdir $ref ${fastqdir}/testsingle.bqB.fastq >${outdir}/.tophat.sbatch
EOF
) > sbatch.tophat

cat sbatch.tophat >tophat.sbatch
sbatch tophat.sbatch

#find . -name 'sbatch.tophat' | xargs -I% sbatch %



#submitted
