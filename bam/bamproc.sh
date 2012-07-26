#Processing bam files
#DEP: rmdups.pl
#DEP: mergebams.pl


###
#rmDups
###
email='alva.rani@scilifelab.se'
projid='b2012046'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
sbatchdir='/proj/b2012046/rani/scripts/gsnap/bamscripts/'
bamfile='/proj/b2012046/rani/analysis/gsnap/'
cd $sbatchdir
find $bamfile -name '*.bam' | awk -F'/' -v execdir=$execdir -v sbatchdir=$sbatchdir -v projid=$projid -v email=$email '{print "perl", execdir"/rmdups.pl", $0, $NF, sbatchdir, projid, email;}' >rmdups.sh
sh rmdups.sh
find . -name '*.sbatch' | xargs -I% sbatch %
#### 

#Check:
find $bamfile -name '*.bam.sorted.nodup' | xargs -I% ls -lrth %
find $bamfile -name '*.rmdups.stderr' | xargs -I% grep -i 'err'
find $bamfile -name '*.rmdups.stderr' | xargs -I% grep -i 'warn'

#rm unsorted bams
find $bamfile -name '*.bam'  | xargs -I% rm %

#rename bamfiles
find $bamfile -name '*.bam.sorted' | xargs -I% echo rename '.bam.sorted' '.sorted.bam' % >cmds.sh
sh cmds.sh
####

###
#Merge and sort bam files
###
projid='b2012046'
email='alva.rani@scilifelab.se'
sbatchdir='/proj/b2012046/rani/scripts/gsnap/mergebams'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
bamfile='/proj/b2012046/rani/analysis/gsnap/'
outdir=$bamfile
time='1:00:00'
cd $sbatchdir
find $bamfile -name '*.sorted.bam' | sed 's/\..*//' | sort -u >samples.list
cat samples.list | xargs -I% basename % | xargs -I% echo perl ${execdir}'/mergebams.pl' % $bamfile $outdir $sbatchdir $projid $email $time >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.sbatch' | xargs -I% sbatch %

#Check:
find ${bamfile}/info/*.stderr | xargs -I% grep -i 'err' %
find ${bamfile}/info/*.mergesams.*.stderr | xargs -I% grep -i 'warn' %
####

###
#Create indexes
###
find $bamfile -maxdepth 1 -name '*.bam' | xargs -I% echo samtools index % > cmd.sh

#Manually add sbatch header to the cmds.sh script
(cat <<EOF
#!/bin/bash -l
#SBATCH -A $projid
#SBATCH -t 30:00
#SBATCH -J indexbam
#SBATCH -p core -n 1
#SBATCH -e ${bamfile}/info/indexbam.jid_%j.stderr
#SBATCH -o ${bamfile}/info/indexbam.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
module load bioinfo-tools
module load samtools/0.1.18
EOF
) >bamindex.sbatchheader
cat bamindex.sbatchheader cmd.sh >bamindex.sbatch
sbatch bamindex.sbatch
