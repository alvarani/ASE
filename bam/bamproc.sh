#Processing bam files
#DEP: rmdups.pl
#DEP: mergebams.pl


###
#rmDups
###
email='alva.rani@scilifelab.se'
projid='b2012046'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
bamdir='/proj/b2012046/edsgard/ase/sim/data/tophat'
sbatchdir='/proj/b2012046/rani/scripts/gsnap/bamscripts'
cd $sbatchdir
find $bamdir -name 'accepted_hits.bam' | awk -F'/' -v execdir=$execdir -v sbatchdir=$sbatchdir -v projid=$projid -v email=$email '{print "perl", execdir"/rmdups.pl", $0, $9, sbatchdir, projid, email;}' >rmdups.sh
sh rmdups.sh
find . -name '*.sbatch' | xargs -I% sbatch %

#Check:
find $bamdir -name '*.bam.sorted.nodup' | xargs -I% ls -lrth %
find $bamdir -name '*.rmdups.stderr' | xargs -I% grep -i 'err'
find $bamdir -name '*.rmdups.stderr' | xargs -I% grep -i 'warn'

#rm unsorted bams
find $bamdir -name 'accepted_hits.bam' | grep -v 'mmseq' | xargs -I% rm %

#rename bamfiles
find $bamdir -name '*.bam.sorted.nodup' | xargs -I% echo rename '.bam.sorted.nodup' '.sorted.nodup.bam' % >cmds.sh
sh cmds.sh


###
#Merge and sort bam files
###
projid='b2012046'
email='alva.rani@scilifelab.se'
sbatchdir='/proj/b2012046/rani/scripts/gsnap/mergebams'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/bam'
bamdir='/proj/b2012046/edsgard/ase/sim/data/tophat'
outdir=$bamdir
time='11:00:00'
cd $sbatchdir
find $bamdir -name '*.sorted.nodup.bam' | sed 's/\..*//' | sort -u >samples.list
cat samples.list | xargs -I% basename % | xargs -I% echo perl ${execdir}'/mergebams.pl' % $bamdir $outdir $sbatchdir $projid $email $time >cmds.sh
sh cmds.sh
find $sbatchdir -name '*.mergesams.sbatch' | xargs -I% sbatch %

#Check:
find ${bamdir}/info/*.mergesams.*.stderr | xargs -I% grep -i 'err' %
find ${bamdir}/info/*.mergesams.*.stderr | xargs -I% grep -i 'warn' %


###
#Create indexes
###
find $bamdir -maxdepth 1 -name '*.bam' | xargs -I% echo samtools index % > cmds.sh

#Manually add sbatch header to the cmds.sh script
(cat <<EOF
#!/bin/bash -l
#SBATCH -A $projid
#SBATCH -t 5:00:00
#SBATCH -J indexbam
#SBATCH -p core -n 1
#SBATCH -e ${bamdir}/info/indexbam.jid_%j.stderr
#SBATCH -o ${bamdir}/info/indexbam.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
module load bioinfo-tools
module load samtools/0.1.18
EOF
) >bamindex.sbatchheader
cat bamindex.sbatchheader cmds.sh >bamindex.sbatch
sbatch bamindex.sbatch
