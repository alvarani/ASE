#Processing bam files
#DEP: rmdups.pl
#DEP: mergebams.pl


###
#rmDups
###
email='alva.rani@scilifelab.se'
projid='b2012046'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
sbatchdir='/proj/b2012046/rani/scripts/gsnap/bamscripts' #/proj/b2012046/rani/scripts/gsnap/umapped --has the scripts for 19 multi mapped reads
bamfile='/proj/b2012046/rani/analysis/gsnap/bamsort/'
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
find  -name '*.bam.sorted' | xargs -I% echo rename '.bam.sorted' '.sorted.bam' % >cmds.sh
sh cmds.sh


#rename bamfiles
find $bamdir -name '*.bam.sorted.nodup' | xargs -I% echo rename '.bam.sorted.nodup' '.sorted.nodup.bam' % >cmds.sh
sh cmds.sh
####

###
#Merge and sort bam files
###
projid='b2012046'
email='alva.rani@scilifelab.se'
sbatchdir='/proj/b2012046/rani/scripts/gsnap/mergebams'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
bamsort='/proj/b2012046/rani/analysis/gsnap/bamsort'
outdir='/proj/b2012046/rani/analysis/gsnap/mergebam'
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
bamsort='/proj/b2012046/rani/analysis/gsnap/mergebam'
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

