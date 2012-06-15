#GSNAP
# I have downloaded the gmap from the homepage and then unziped and extracted
wget 'http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-05-15.tar.gz'
gunzip gmap-gsnap-2012-06-12.tar.gz
tar xvf./* gmap-gsnap-2012-06-12.tar
tar -xvf gmap-gsnap-2012-06-12.tar
#Made some changes in the config file
#Specified the paths
./configure
make 
make installg
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
cat Ailuropoda_melanoleuca.ailMel1.67.gtf | /bubo/home/h24/alvaj/gsnap-2012-06-12/util/gtf_splicesites > at.splicesites
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
# wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/snp130.txt.gz
# unzipped by,
gunzip -c snp130.txt.gz | ./dbsnp_iit -w 3 -e snp130Exceptions.txt.gz > snpfile.txt
## changing the snpfile to .iit format file by,
cat snpfile.txt | /bubo/home/h24/alvaj/gmap-2012-05-15/src/iit_store -o snpfile
## after cpying the itt file into hg19.maps folder , and the programm to find the indictaed file from the .maps folder  by the folowwing command(used at1.iit file)
/bubo/home/h24/alvaj/gmap-2012-05-15/src/snpindex -d hg19  -v at1

# for intron sites
cat Homo_sapiens.GRCh37.59.gtf |/bubo/home/h24/alvaj/glob/annotation/gsnap/gmap-2012-05-15/util/gtf_introns > snp.introns

(cat <<EOF
#!/bin/bash -l
#SBATCH -A $projid
#SBATCH -t 55:00
#SBATCH -J haploref
#SBATCH -p core -n 1
#SBATCH -e $outdir/phasing/log/haploref.jid_%j.stderr
#SBATCH -o $outdir/phasing/log/haploref.jid_%j.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=$email
export PATH=$PATH:/bubo/home/h26/edsgard/opt/mmseq_0.10.0
cd ${outdir}/phasing
haploref.rb ${BOWTIE_INDEXES}$TRANSCRIPT_FASTA ${BOWTIE_INDEXES}$GFF_FILE SNPinfo sample.hap > sample.fa
EOF
) >sbatch.template
cat sample.list | xargs -I% echo cat sbatch.template "| sed 's/sample/"%"/g' >" %.haploref.sbatch >cmds.sh
sh cmds.sh
find . -name '*.haploref.sbatch' | xargs -I% sbatch %

##****Kalkyl****
# assigning the job for some other interactive level
interactive -p devel -t 1:00:00 -A b2012046
# to see om whic node I am runing
srun hostname -s |sort -u
# to create same directory for all the nodes which I am working with
 srun -N 4 -n 4 mkdir /scratch/$SLURM_JOB_ID/indata
# Copy indata for my_program to the local directories
 srun -N 4 -n 4 cp -r ~/glob/indata/* /scratch/$SLURM_JOB_ID/indata
# to see the queue line of the job 
squeue -u alvaj


# to rename the file 
cat snp.splicesite | awk '{print $1, "chr"$2, $3, $4;}' >snp.splicesite.uscschr
# tried for a short seqeunce from the fastq
/bubo/home/h24/alvaj/gmap-2012-05-15/src/gsnap -d hg19 ./fastq >shortfastaseq
# to copy  a .sh file from the local computer to the cloud...
scp g.sh alvaj@kalkyl.uppmax.uu.se:/bubo/home/h24/alvaj/glob/annotation/gsnap/g.sh

/bubo/home/h24/alvaj/gmap-2012-05-15/src/gsnap -d hg19 -s./splicesitesfile /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/1_unstim.readp_1.2.ill.fastq.1.filter.fastq >output1

# for paired end rna seq reads
/bubo/home/h24/alvaj/gmap-2012-05-15/src/gsnap -d hg19 -s./splicesitesfile /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/1_unstim.readp_1.2.ill.fastq.1.filter.fastq /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/ 9_unstim.readp_1.2.ill.fastq.1.filter.fastq >output2


(nohup prg args >my.out) >& my.err &
sbatch .script
srun -A b21010 -t 1:00:00 -p devel prg args >my.out

##
#RAN for small seq and is working fine!!!
annotdir='/bubo/home/h24/alvaj/glob/annotation/gsnap'
outdir='/proj/b2012046/rani/analysis/gsnap'
datadir='/proj/b2012046/rani/data/fastq'
splicefile=${annotdir}/'splicesiteschro'
fastqfile1=${datadir}/test.n100.rp1.fastq
fastqfile2=${datadir}/test.n100.rp2.fastq
refdir=${annotdir}/gmapdb
ref=hg19
snpfile=dbsnp135
snpdir=${annotdir}/dbsnp
samdir='/bubo/home/h24/alvaj/glob/annotation/gsnap/samtools/samtools-0.1.6'

cd $outdir
export PATH=$PATH:/bubo/home/h24/alvaj/opt0/gmap-2012-06-12/bin

(gsnap -D $refdir -d $ref -A sam -s $splicefile -V $snpdir -v $snpfile --quality-protocol=illumina $fastqfile1 $fastqfile2 | ${samdir}/samtools view -Sb - > ${outdir}/test.n100.bam12) >& test.err12 &

| ${samdir}/samtools view -Sb test.sam > test.sam
| ${samdir}/samtools view -S sam -b bam
samtools view -bT hg.fa -o xxx.bam xxx.sam
samdir='/bubo/home/h24/alvaj/glob/annotation/gsnap/samtools/samtools-0.1.6'

## github commit// to push the changes from my computer to github
#when a fatal error occurs such as index lock already running
 rm -f ./.git/index.lock

# github to add file from computer towards repo
git add run.tophat.sh 
 git commit -m 'The edited tophat.sh file' run.tophat.sh
git remote add origin https:https://github.com/alvarani/ASE.git
git push origin master

#kill all slurm jobs--if you want to specify the number of jobs eg: 'NR >32 {print $1;}'
squeue -u alvaj | grep -v JOB | awk '{print $1;}' | xargs -I% scancel %


# stream editor SED an example:to print the lines containg s/1.filter.fastq// and 2.filter.fastq'....
cat fastq.files | sed 's/1.filter.fastq//' | grep -v '2.filter.fastq' >fastqfiles.prefix


# to delete a file with only a particular format within a directory,here to delete the files with .sam format
find . -type f -name "*.sam" -exec rm -f {} \;



# it helps to change the suffix from fastqfiles10.prefix to sample1 and sample2  
cat fastqfiles10.prefix | xargs -I% basename % | xargs -I% echo cat sbatch10.template "| sed 's/sample/"%"/g' >" %gsnap.sbatch >cmds10.sh



# to count the number of files in a directory, where targetdir is the name of directory
ls -1 targetdir | wc -l
# which is more advanced
find Downloads -maxdepth 1  -type f | wc -l



# listing the jobs in human readable form
ls -ltrh
#to convet from sam to bam
samdir= '/bubo/home/h24/alvaj/glob/annotation/gsnap/samtools/samtools-0.1.6/samtools'
outdir='/proj/b2012046/rani/analysis/gsnap'
find . -type f -name "*.sam" | ${samdir}samtools view -bS  > ${outdir}/samplebam


#example sed 
sed 's/day/night/' <old >new

s	  Substitute command
/../../	  Delimiter
day	  Regular Expression Pattern Search Pattern
night	  Replacement string , ie replace day with night
# anothe rexample with echno
echo Sunday | sed 's/day/night/'

# the codes for bam file is in this one ASE
/bubo/home/h24/alvaj/glob/code/ASE

