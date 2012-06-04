#Stats on bam files
#DEP: coverage.pl
#DEP: bedtoolshist2pdf.R


###
#Calculate number of mapped reads (More than the input number of reads, due to non-unique mappings!)
###
module load bioinfo-tools
module load samtools/0.1.16
find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -name 'accepted_hits.bam.sorted.nodup' -exec sh -c "samtools flagstat {} >{}.samtools.flagstat.out" \; &
find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -name '*.samtools.flagstat.out' | xargs -I% awk 'NR == 1 {print FILENAME, $1;}' % >nmapped.reads


###
#Calculate coverage
###

#BEDtools (see coverage.pl)
prg='/bubo/home/h26/edsgard/glob/code/ngs_pipes/hs_pipe/coverage.pl'
find '/proj/b2011075/analysis/hs_pipe/outdata/run1' -name 'merged.bam' | awk -F'/' -v ods="$ods" -v projid="$projid" -v em="$email" -v prg="$prg" '{print "perl", prg, $0, $8, ods, projid, em;}' >bedtoolscov.sh


###
#Generate pdf of the coverage
###
find . -name 'ccds.bedtools.out2' -exec sh -c 'grep ^all {} >{}.all' \;
find . -name 'ensembl.bedtools.out' -exec sh -c 'grep ^all {} >{}.all' \;
find . -name 'genome.bedtools.out2' -exec sh -c 'grep ^genome {} >{}.all' \;
find . -name 'ccds.*2.all' | zip ccds.bedtools.all.zip -@
find . -name 'genome.*2.all' | zip genome.bedtools.all.zip -@
find . -name 'ccds*2.all' | awk -F'/' '{print "mv", $0, "coverage/"$2"."$4;}' >mv.sh
find . -name 'genome*2.all' | awk -F'/' '{print "mv", $0, "coverage/"$2"."$4;}' >mv.sh
find '/Users/danieledsgard/Dropbox/postdoc/projects/ase/res/hs_pipe/run1' -name '*ccds.*2.all' | sed 's/.*\///' | awk -F'.' '{print $1;}' >samples
find '/Users/danieledsgard/Dropbox/postdoc/projects/ase/res/hs_pipe/run1' -name '*ccds.*2.all' >files
paste samples files >ccds.histfiles
find '/Users/danieledsgard/Dropbox/postdoc/projects/ase/res/hs_pipe/run1' -name '*genome.*2.all' | sed 's/.*\///' | awk -F'.' '{print $1;}' >samples
find '/Users/danieledsgard/Dropbox/postdoc/projects/ase/res/hs_pipe/run1' -name '*genome.*2.all' >files
paste samples files >genome.histfiles

bedtoolshist2pdf.R ccds.histfiles

#NB: no pair after read-QC, rmdups, but not: unique, MQ, properly paired. What is M (mismatches)set to?

