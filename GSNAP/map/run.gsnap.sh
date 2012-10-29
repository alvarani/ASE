
# download the software and DB
wget 'http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-06-12.tar.gz'
wget 'http://research-pub.gene.com/gmap/genomes/hg19.tar'

gsnapexecdir='/bubo/home/h24/alvaj/opt0/gmap-2012-06-12'

#SNPs
###
# downloaded the needed file snp135common, 
snpdir='/bubo/home/h24/alvaj/glob/annotation/gsnap'
gsnapexecdir='/bubo/home/h24/alvaj/opt0/gmap-2012-06-12'
snpfile=dbsnp135

cd $snpdir
## done again saved in annotation/gsnap for latest version of gsanp

wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp135Common.txt.gz'



#reformat the dbsnp 
gunzip -c snp135Common.txt.gz | ${gsnapexecdir}/dbsnp_iit  -w 3  > ${snpfile}.txt &
gunzip -c snp135Common.txt.gz | ${gsnapexecdir}/dbsnp_iit  -w 3  > dbsnp135.txt &


(cat ${snpfile}.txt |${gsnapexecdir}/iit_store -o $snpfile >${snpfile}.iitstore.out) >& ${snpfile}.iitstore.err &

(snpindex -D $refdir -d $ref -V $snpdir -v $snpfile ${snpfile}.iit >snpindex.out) >& snpindex.err &


# Splice sites , to generate a splice site index

cat Homo_sapiens.GRCh37.59.gtf |${gsnapexecdir}/gtf_splicesites > snp.splicesiteschr


# renamed file is splicesitechr
cat snp.splicesiteschr | awk '{print $1, "chr"$2, $3, $4;}' >splicesiteschr 




## processing it to map file (changing to .iit file)
cat splicesiteschr | ${gsnapexecdir}/iit_store -o splicesiteschro


###
## Run gsnap
###

#set path for for each of these variables 
gsnapexecdir='/bubo/home/h24/alvaj/opt0/gmap-2012-06-12'
annotdir='/bubo/home/h24/alvaj/glob/annotation/gsnap'
scriptdir='/proj/b2012046/rani/scripts/gsnap/mapscripts'
fastqdir='/proj/b2012046/rani/data/synt/fastqfilt/oneGchunks'
outdir='/proj/b2012046/rani/analysis/gsnap'
splicefile=${annotdir}/'splicesiteschro'
refdir=${annotdir}/gmapdb
ref=hg19
snpdir=${annotdir}/dbsnp
snpfile=dbsnp135
projid='b2012046'
email='alva.rani@scilifelab.se'

#Rename fastq files and creating a link to the fastq  data directory
cd $fastqdir
ln -s /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/oneGchunks/*.fastq .
#append readp_1 as suffix
ls -1 *.readp_1.*.fastq | xargs -I% mv % %.readp_1
#rm first occurence of readp_1
rename .readp_1.part .part *.readp_1
#append readp_2 as suffix
ls -1 *.readp_2.*.fastq | xargs -I% mv % %.readp_2
#rm first occurence of readp_2
rename .readp_2.part .part *.readp_2


# creation of the batch scripts for the paired end reads
# choosing the prefix of the samples
cd $scriptdir
find $fastqdir -maxdepth 1 -name '*fastq.readp_1' | sed 's/.readp_1//' >fastq.files.prefix

samples=(`cat fastq.files.prefix`)

for fastqfile in ${samples[@]}
do
sbatchfile=`basename $fastqfile`
echo $sbatchfile
(cat <<EOF
#!/bin/bash -l
#SBATCH -A $projid
#SBATCH -t 30:00:00
#SBATCH -J gsnap
#SBATCH -p node -n 1
#SBATCH -e $outdir/log/gsnap.${sbatchfile}.jid_%j.stderr
#SBATCH -o $outdir/log/gsnap.${sbatchfile}.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
export PATH=$PATH:${gsnapexecdir}/bin
cd ${fastqdir}
module load bioinfo-tools
module load samtools/0.1.18
 gsnap -D $refdir -d $ref  -A sam -s $splicefile -V $snpdir -v $snpfile --quality-protocol=illumina  --npaths=1  -Q --quiet-if-excessive ${fastqfile}.readp_1 ${fastqfile}.readp_2 | samtools view -Sb - >${outdir}/${sbatchfile}.bam
EOF
)>${sbatchfile}.sbatch
done
find . -name '*.fastq.sbatch' | xargs -I% sbatch %

