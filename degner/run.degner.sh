
#Simulate data a la Degner
#DEP: get.overlaps.v3.perl (from Degner)


###
#Format input variants
###
#id chr pos strand ref alt
#EX: 10.100000020 10 100000020 + A G

vcfdir='/proj/b2012046/edsgard/ase/sim/data/varcalls'
vars=${vcfdir}/vars2degner.dp10.txt
cat ${vcfdir}/*.hetvars.alleles | sort -u | awk -F'\t' '{print $1"."$2, $1, $2, "+", $3, $4;}' >$vars

#Check:
wc -l vars2degner.dp10.txt #112901


###
#Create simulated read set
###
#TBD: Add noise.
#TBD: NB: Disregards if >1 SNP covering a read!
#TBD: PE reads not simulated, only SEs!
degneroutdir='/proj/b2012046/edsgard/ase/sim/data/degner'
fastqoutfile=${degneroutdir}/mpileup.dp10.len100.bqB.fastq
vcfdir='/proj/b2012046/edsgard/ase/sim/data/varcalls'
vars=${vcfdir}/vars2degner.dp10.txt
projid='b2010035'
time='4:57:00'
logdir=${degneroutdir}/log
exedir='/bubo/home/h26/edsgard/glob/code/ase/degner'
refdir='/bubo/home/h26/edsgard/glob/annotation/bowtie'
readlen='100'
readqual='BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'

mkdir -p $logdir
srun -A $projid -p node -t $time -e ${logdir}/sim.err -o $fastqoutfile perl ${exedir}/get.overlaps.v3.perl $vars $refdir $readlen $readqual &
#Status: submitted


###
#Map with Tophat
###
#See map/run.tophat.sh
#You'll need to amend it such that it takes single end reads instead of paired end reads!


###
#rmDups
###
#Skip rmdups


###
#Merge and sort bam files
###
#TBD: Not sure that you need this. But in that case see bam/bamproc.sh


###
#MPILEUP (no varcall)
###
#See basecount/basecount.sh


###
#Mapping bias modification a la Montgomery
###
#I'll get back to you with code for that once you've reached this point:)
