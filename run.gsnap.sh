#! /bin/bash -l
# download the software and DB
wget 'http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-06-06.tar.gz'
wget 'http://research-pub.gene.com/gmap/genomes/hg19.tar'

gsnapexecdir='/bubo/home/h24/alvaj/opt0/gmap-2012-06-06'

#SNPs
###
# downloaded the needed file snp134,snp135common, 
cd $snpdir
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp135Common.txt.gz'
## done again saved in annotation/gsnap for latest version of gsanp
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp135.txt.gz'

#reformat---from here I have to run for the new version
gunzip -c snp135Common.txt.gz | ${gsnapexecdir}/dbsnp_iit  -w 3  > ${snpfile}.txt &
gunzip -c snp135Common.txt.gz | ${gsnapexecdir}/dbsnp_iit  -w 3  > dbsnp135.txt &
## done and save in annotation/gsnap |I used it here for new version)

(cat ${snpfile}.txt |${gsnapexecdir}/iit_store -o $snpfile >${snpfile}.iitstore.out) >& ${snpfile}.iitstore.err &
##done for new version

(snpindex -D $refdir -d $ref -V $snpdir -v $snpfile ${snpfile}.iit >snpindex.out) >& snpindex.err &
/bubo/home/h24/alvaj/opt0/gmap-2012-06-06/src/snpindex

###
# Splice sites , to generate a splice site index---***splicesite***
###
cat Homo_sapiens.GRCh37.59.gtf |${gsnapexecdir}/gtf_splicesites > snp.splicesiteschr
#(---done)

cat snp.splicesiteschr | awk '{print $1, "chr"$2, $3, $4;}' >splicesiteschr 
#(has to do it agagin)******IMP****since rewritten over it

# renamed file is splicesitechr

## processing it to map file (changing to .iit file)---done for snp.splicesitechr
cat splicesiteschr | ${gsnapexecdir}/iit_store -o splicesiteschro
## done

###
## Fastq files
###

#set path for for each of these variables 
gsnapexecdir='/bubo/home/h24/alvaj/opt0/gmap-2012-06-06'
annotdir='/bubo/home/h24/alvaj/glob/annotation/gsnap'
scriptdir='/proj/b2012046/rani/scripts/gsnap'
fastqdir='/proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt'
outdir='/proj/b2012046/rani/analysis/gsnap'

splicefile=${annotdir}/'splicesiteschro'
refdir=${annotdir}/gmapdb
ref=hg19
snpdir=${annotdir}/dbsnp
snpfile=dbsnp135
projid='b2012046'
email='alva.rani@scilifelab.se'


# all files apart from those with '.S' extension
cd $scriptdir
find $fastqdir -maxdepth 1 -name '*fastq' | grep -v '\.S\.' >fastq.files

#Create sbatch scripts
cat fastq.files | sed 's/1.filter.fastq//' | grep -v '2.filter.fastq' >fastqfiles.prefix

(cat <<EOF
#!/bin/bash -l
#SBATCH -A $projid
#SBATCH -t 35:00:00
#SBATCH -J gsnap
#SBATCH -p node -n 1
#SBATCH -e $outdir/log/gsnap.samplejid_%j.stderr
#SBATCH -o $outdir/log/gsnap.samplejid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
export PATH=$PATH:${gsnapexecdir}/bin
cd ${fastqdir}
gsnap -D $refdir -d $ref -A sam -s $splicefile -V $snpdir -v $snpfile --quality-protocol=illumina sample1.filter.fastq sample2.filter.fastq >${outdir}/samplesam
EOF
) >sbatch.template
cat fastqfiles.prefix | xargs -I% basename % | xargs -I% echo cat sbatch.template "| sed 's/sample/"%"/g' >" %gsnap.sbatch >cmds.sh
sh cmds.sh
find . -name '*.gsnap.sbatch' | xargs -I% sbatch %
#status: submitted




