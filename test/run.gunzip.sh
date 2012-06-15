#GSNAP
# I have downloaded the gmap from the homepage and then unziped and extracted
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-05-15.tar.gz
gunzip gmap-gsnap-2012-05-15.tar.gz
tar xvf./* gmap-gsnap-2012-05-15.tar
#Made some changes in the config file
#Specified the paths
./configure
make 
make install
# I have downloaded the database and extracted it
wget http://research-pub.gene.com/gmap/genomes/hg19.tar
tar xvf  hg19.tar
##optional if need to run with older versions
 ln -s hg19.ref153positions hg19.ref12153.positions
#the  fasta file used here
wget 'http://bgx.org.uk/software/Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa.gz'
gunzip Homo_sapiens.GRCh37.64.dna.toplevel.ref.fa.gz

#How to process a genome for GMAP...Started with an example file 'g'
# Ran an example ,downloaded a fasta file and saved as g in utils
#Ran coords.txt
./fa_coords g
#then Ran gmap_setup
./gmap_setup -d hg19 g
# tried the file from gmapdb
./md_coords /bubo/home/h24/alvaj/gmap-2012-05-15/gmapdb/hg19/hg19.contig


