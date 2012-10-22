
###Test for 19 sam file
#created a errorlist of 19 files
# taking the the reads appreaed more than twice and taking it in a order 
cat errorfiles.list | xargs -I% echo "samtools view "%" | awk '{print \$1;}' | sort | uniq -c >"%.reads2nmapped >cmds.sh



errordir='/proj/b2012046/rani/analysis/gsnap/bamsort/errortest'

ls -1 ${errordir}/*.reads2nmapped >reads2nmapped.files

cat reads2nmapped.files | xargs -I% basename % | xargs -I% echo "awk '\$1 == 2 {print \$0;}'" ${errordir}/% " | wc -l >" %.n.uniq >cmds.sh
cat reads2nmapped.files | xargs -I% basename % | xargs -I% echo "awk '\$1 != 2 {print \$0;}'" ${errordir}/% " | wc -l >" %.n.nonuniq >cmds2.sh


bamfiles='/proj/b2012046/rani/analysis/gsnap/bamsort/errortest/errorbamfiles.list'
errordir='/proj/b2012046/rani/analysis/gsnap/bamsort/errortest'

#Extract multimapped reads
ls -1 ${errordir}/*.reads2nmapped >reads2nmapped.files
cat reads2nmapped.files | xargs -I% basename % | xargs -I% echo "awk '\$1 != 2 {print \$0;}'" ${errordir}/% "| awk '{print \$2;}'>" %.nonuniq >cmds.sh
sh cmds.sh &

#TBD: Convert all bams to sams
bamdir='/proj/b2012046/rani/analysis/gsnap/bamindex'
sam=${bamdir}/2_LPS.part_2ab.sam

#Rm multimapped reads
ls -1 *d.nonuniq | awk -F'.' -v OFS='.' '{print $1, $2;}' >errorfiles.prefix.list
cat errorfiles.prefix.list | xargs -I% echo "perl rm.multimappedreads.pl " %.sam %.fastq.bam.reads2nmapped.nonuniq ">"%.uniqmapped.sam >cmds.sh

## How to add a word or letter to the first column of the file
 cat CCDS-hg19.chrstripped.bed | awk -Fchr '{print "chr"$1"";}' >CCDS-hg19.chrstripped.Nbed
## How to remove lines from 1-24 in afile
cat hg19.chrsizes.chrstripped | sed '1,24d' >body.chrstripped

