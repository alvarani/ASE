Stats on bam files
#DEP: coverage.pl
#DEP: bedtoolshist2pdf.R


###
#Calculate number of mapped reads (More than the input number of reads, due to non-unique mappings!)
###
module load bioinfo-tools
module load samtools/0.1.16
find '/proj/b2012046/rani/analysis/gsnap/bamsort' -name '*.sorted.nodup.bam' -exec sh -c "samtools flagstat {} >{}.samtools.flagstat.out" \; &
find '/proj/b2012046/rani/analysis/gsnap/bamsort' -name '*.samtools.flagstat.out' | xargs -I% awk 'NR == 1 {print FILENAME, $1;}' % >nmapped.reads


###
#Calculate coverage
###
projid='b2012046'
email='alva.rani@scilifelab.se'
#BEDtools (see coverage.pl)
prg='/bubo/home/h24/alvaj/glob/code/ASE/bam/coverage.pl'
find '/proj/b2012046/rani/analysis/gsnap/mergebam' -name '*.merged.bam' | awk -F'/' -v ods="$ods" -v projid="$projid" -v em="$email" -v prg="$prg" '{print "perl", prg, $0, $8,ods, projid, em;}' >bedtoolscov.sh
sh bedtoolscov.sh

find . -name '*.sbatch2' | xargs -I% sbatch %

###
#Generate pdf of the coverage for the ccds file:
###
find . -name 'ccds.bedtools.out2' -exec sh -c 'grep ^all {} >{}.all' \;
find . -name 'ccds.*2.all' | sed 's/\.*\///' | awk -F'.' '{print $1;}' >samples
find /proj/b2012046/rani/analysis/degner -name '*ccds.*2.all' >files
paste samples files >ccds.histfiles
find . -name 'ccds*2.all' | awk -F'/' '{print "mv", $0, "coverage/"$2"."$4;}' >mv.sh

bedtoolshist2pdf.R ccds.histfiles

#NB: no pair after read-QC, rmdups, but not: unique, MQ, properly paired. What is M (mismatches)set to?


