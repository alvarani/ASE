#Processing bam files
#DEP: rmdups.pl
#DEP: mergebams.pl


###
#rmDups
###
email='alva.rani@scilifelab.se'
projid='b2012046'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
sbatchdir='/proj/b2012046/rani/scripts/degner/bamscripts'
bamfile='/proj/b2012046/rani/analysis/degner'
cd $sbatchdir
find $bamfile -name 'accepted_hits.bam' | awk -F'/' -v execdir=$execdir -v sbatchdir=$sbatchdir -v projid=$projid -v email=$email '{print "perl", execdir"/rmdups.pl", $0, $7, sbatchdir, projid, email;}' >rmdups.sh
sh rmdups.sh
find . -name '*.sbatch' | xargs -I% sbatch %
####
find . -name '*.sorted' | wc -l

#Check:
find $bamfile -name '*.bam.sorted.nodup' | xargs -I% ls -lrth %
find $bamfile -name '*.rmdups.stderr' | xargs -I% grep -i 'err'
find $bamfile -name '*.rmdups.stderr' | xargs -I% grep -i 'warn'

#rm unsorted bams
find $bamfile -name '*.bam'  | xargs -I% rm %

#rename bamfiles
find  -name '*.bam.sorted' | xargs -I% echo rename '.bam.sorted' '.sorted.bam' % >cmds.sh
sh cmds.sh

#rm unsorted bams
find $bamdir -name 'accepted_hits.bam' | grep -v 'mmseq' | xargs -I% rm %

#rename bamfiles
find $bamdir -name '*.bam.sorted.nodup' | xargs -I% echo rename '.bam.sorted.nodup' '.sorted.nodup.bam' % >cmds.sh
sh cmds.sh
####

###
#Merge and sort bam files
###
projid='b2012046'
email='alva.rani@scilifelab.se'
sbatchdir='/proj/b2012046/rani/scripts/degner/mergebamScripts'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
bamsort='/proj/b2012046/rani/analysis/degner/'
outdir='/proj/b2012046/rani/analysis/degner/mergebam'
time='11:00:00'
cd $sbatchdir
find $bamsort -name '*.sorted.nodup.bam' | sed 's/\..*//' | sort -u >samples.list
cat samples.list | xargs -I% basename % | xargs -I% echo perl ${execdir}'/mergebams.pl' % $bamsort $outdir $sbatchdir $projid $email $time >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.mergesams.sbatch' | xargs -I% sbatch %



## Submitted
#Check:
find ${bamfile}/info/*.stderr | xargs -I% grep -i 'err' %
find ${bamfile}/info/*.mergesams.*.stderr | xargs -I% grep -i 'warn' %
####

###
#Create indexes
###
bamsort='/proj/b2012046/rani/analysis/degner/mergebam'
find $bamsort -maxdepth 1 -name '*.bam' | xargs -I% echo samtools index % > cmd.sh
#Manually add sbatch header to the cmds.sh script
(cat <<EOF
#!/bin/bash -l
#SBATCH -A $projid
#SBATCH -t 3:00:00
#SBATCH -J indexbam
#SBATCH -p core -n 1
#SBATCH -e ${bamsort}/info/indexbam.jid_%j.stderr
#SBATCH -o ${bamsort}/info/indexbam.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
module load bioinfo-tools
module load samtools/0.1.18
EOF
) >bamindex.sbatchheader
cat bamindex.sbatchheader cmd.sh >bamindex.sbatch
sbatch bamindex.sbatch

############################
## Done all steps here without errors

