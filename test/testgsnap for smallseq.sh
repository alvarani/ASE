


#RUN!for small test  seq

annotdir='/bubo/home/h24/alvaj/glob/annotation/gsnap'
splicefile=${annotdir}/'splicesiteschro'
snpfile=dbsnp135
snpdir=${annotdir}/dbsnp
outdir='/proj/b2012046/rani/analysis/gsnap'
datadir='/proj/b2012046/rani/data/fastq'
fastqfile1=${datadir}/test.n100.rp1.fastq
fastqfile2=${datadir}/test.n100.rp2.fastq
refdir=${annotdir}/gmapdb
ref=hg19

cd $outdir
export PATH=$PATH:/bubo/home/h24/alvaj/opt/gmap-2012-06-02/bin
(gsnap -D $refdir -d $ref -A sam -s $splicefile -V $snpdir -v $snpfile --quality-protocol=illumina $fastqfile1 $fastqfile2 >${outdir}/test.n100.sam) >& test.err &

### worked and given output

