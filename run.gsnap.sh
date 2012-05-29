### *****-----Started with SNP tolerant alignment----****

wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp135Common.txt.gz'
gunzip snp135Common.txt.gz

#changing the format from .txt to iit 
cat snpfile.txt | /bubo/home/h24/alvaj/gmap-2012-05-15/src/iit_store -o snpfile
# saved in hg19.maps folder

#----Detecting known and novel splice sites in GSNAP----

## download the fastaq file from the directory
fastqdir='/proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt'

# to download all files else with '.S' extension
find $fastqdir -maxdepth 1 -name '*fastq' | grep -v 'S' >fastq.files

# Splice sites , to generate a splice site index
cat Homo_sapiens.GRCh37.59.gtf |/bubo/home/h24/alvaj/glob/annotation/gsnap/gmap-2012-05-15/util/gtf_splicesites > snp.splicesites

## processing it to map file (changing to .iit file)
cat snp.splicesites |/bubo/home/h24/alvaj/glob/annotation/gsnap/gmap-2012-05-15/src/iit_store -o splicesitesfile

/bubo/home/h24/alvaj/gmap-2012-05-15/src/gsnap -d hg19 -s./splicesitesfile /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/1_unstim.readp_1.2.ill.fastq.1.filter.fastq >output1

