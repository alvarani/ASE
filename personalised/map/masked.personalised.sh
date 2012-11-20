

#!/bin/bash -l
#SBATCH -A b2012046
#SBATCH -t 5:00:00
#SBATCH -J personlised
#SBATCH -p node -n 8
#SBATCH -e /bubo/home/h24/alvaj/glob/code/ASE/personalised/1_LPS.readp_1.1.ill.fastq.err
#SBATCH -o /bubo/home/h24/alvaj/glob/code/ASE/personalised/1_LPS.readp_1.1.ill.fastq.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=alva.rani@scilifelab.se
module load bioinfo-tools
module load samtools/0.1.18
module load tophat/2.0.4
tophat --solexa1.3-quals -p 8 -o /proj/b2012046/rani/analysis/degner/personalised/1_LPS.readp_1.1.ill.fastq --GTF /bubo/home/h26/edsgard/glob/annotation/human/Homo_sapiens.GRCh37.59.gtf -r '-36' --mate-std-dev 75 -g 1  -x 1  -M  /bubo/home/h24/alvaj/glob/annotation/degner/personalised/Homo_sapiens.maskedRefe.personalised.GRCh37.57.fa /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/1_LPS.readp_1.1.ill.fastq.1.filter.fastq /proj/b2012046/edsgard/ase/sim/data/synt/fastqfilt/1_LPS.readp_1.1.ill.fastq.2.filter.fastq


#Bowtie bulid
/bubo/home/h24/alvaj/glob/annotation/degner/reference/Homo_sapiens.maskedRefe.personalised.GRCh37.57.fa

#!/bin/bash -l
#SBATCH -A b2012046
#SBATCH -t 5:00:00
#SBATCH -J pwersonlised
#SBATCH -p node -n 8
#SBATCH -e /bubo/home/h24/alvaj/glob/code/ASE/personalised/bulid.err
#SBATCH -o /bubo/home/h24/alvaj/glob/code/ASE/personalised/bulid.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=alva.rani@scilifelab.se
module load bioinfo-tools
module load bowtie/0.12.6
bowtie-build -f /bubo/home/h24/alvaj/glob/annotation/degner/personalised/Homo_sapiens.maskedRefe.personalised.GRCh37.57.fa /bubo/home/h24/alvaj/glob/annotation/degner/reference/

