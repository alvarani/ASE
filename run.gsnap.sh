#GSNAP
# I have downloaded the gmap from the homepage and then unziped and extracted
wget 'http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-05-15.tar.gz'
gunzip gmap-gsnap-2012-05-15.tar.gz
tar xvf./* gmap-gsnap-2012-05-15.tar
#Made some changes in the config file
#Specified the paths
./configure
make 
make install
# I have downloaded the database and extracted it
wget 'http://research-pub.gene.com/gmap/genomes/hg19.tar'
tar xvf  hg19.tar
##optional if need to run with older versions
 ln -s hg19.ref153positions hg19.ref12153.positions
#the  fasta file used here
wget 'http://bgx.org.uk/software/Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa.gz'
gunzip Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa.gz

wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp135.txt.gz'

#How to process a genome for GMAP...Started with an example file 'g'
# Ran an example ,downloaded a fasta file and saved as g in utils
#Ran coords.txt
./fa_coords g
#then Ran gmap_setup
./gmap_setup -d hg19 g
# tried the file from gmapdb
./md_coords /bubo/home/h24/alvaj/gmap-2012-05-15/gmapdb/hg19/hg19.contig

#Building map files
# here I used an example file 'Ailuropoda_melanoleuca.ailMel1.67.gtf',first I concatinated this file with a utility program gtf_splicesites by his command and stored in a folder at.splicesites
cat Ailuropoda_melanoleuca.ailMel1.67.gtf | /bubo/home/h24/alvaj/gmap-2012-05-15/util/gtf_splicesites > at.splicesites
## then rename at.splicesites to at1, the folowwing comand
cat at.splicesites |./iit_store -o at1

## Runing GSANP, the following command helps to run the fastaq file towards the genome hg19
./gsnap -d hg19 x

# GSNAP gives the difference between site level and intron level splicing based on the format of input file, to perform site level splicing the file should be in a special format which can be created by the folowwing commnad
 cat Ailuropoda_melanoleuca.ailMel1.67.gtf |/bubo/home/h24/alvaj/gmap-2012-05-15/util/gtf_splicesites >foo.splicesites
# for intron level the format can be created by,
cat Ailuropoda_melanoleuca.ailMel1.67.gtf |/bubo/home/h24/alvaj/gmap-2012-05-15/util/gtf_introns >foo.introns
## once it is created you can process for map files by the following command,
cat foo.splicesites | ./iit_store -o splicesitefile
## Example for SNP tolerant alignment in GSNAP
 wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/snp130.txt.gz
# unzipped by,
gunzip -c snp130.txt.gz | ./dbsnp_iit -w 3 -e snp130Exceptions.txt.gz > snpfile.txt
## changing the snpfile to .iit format file by,
cat snpfile.txt | /bubo/home/h24/alvaj/gmap-2012-05-15/src/iit_store -o snpfile
## after cpying the itt file into hg19.maps folder , and the programm to find the indictaed file from the .maps folder  by the folowwing command(used at1.iit file)
/bubo/home/h24/alvaj/gmap-2012-05-15/src/snpindex -d hg19  -v at1


### -----Started with SNP tolerant alignment----
 wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp135Common.txt.gz'
gunzip snp135Common.txt.gz


