#Map with tophat
#DEP: tophat.gensbatch.pl


###
#Tophat
###
#Longest run: 5h38min, 6G file for each readpair: 8_LPS.readp_1.2.ill.fastq
email='alva.rani@scilifelab.se'
projid='b2012046'
execdir='/bubo/home/h24/alvaj/glob/code/ASE/map'
fastqdir='/proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt'
fastqsuffix='.filter.fastq'
outdir='/proj/b2012046/rani/analysis/degner'
sbatchdir='/proj/b2012046/rani/scripts/degner'

threads='8'
ref='/bubo/home/h24/alvaj/glob/code/ASE/degner/bowtie-0.12.8/indexes/ref/Ref.fa '

time='23:00:00'

threads='8'
annotdir='/bubo/home/h26/edsgard/glob/annotation/human/'
gtf=${annotdir}'Homo_sapiens.GRCh37.59.gtf'
time='23:00:00' #Longest run: 7h10min: Input: ~7.3G per readpair fastq file
readlen='100'
isizefile='/proj/b2011075/data/rnaseq/bioanalyzer/isizes.only.tab'
isizedev='75' #very rough approx from the bioanalyzer plots. Also, we use trimming which further adds to the deviation.

cd $fastqdir
## creating a sample
ls *.fastq | awk '{split($0,b,"."); print b[1]}' > ${sbatchdir}/samples.list
cd $sbatchdir
rm cmds.sh

samples=(`cat samples.list`)
for sample in ${samples[@]}
do
find $fastqdir -maxdepth 1 -name ${sample}'*' | sed 's/[12].filter.fastq//' | grep -v '.S.filter.fastq' | sort -u >fastqfiles.prefix.list
    cat fastqfiles.prefix.list | xargs -I% echo perl ${execdir}/'tophat.gensbatch.pl' $sample % $fastqsuffix $outdir $threads $gtf $ref $isizefile $sbatchdir $projid $email $time $readlen $isizedev >>cmds.sh
done
sh cmds.sh
find . -name '*.tophat.sbatch' | xargs -I% sbatch %

## Submitted again


