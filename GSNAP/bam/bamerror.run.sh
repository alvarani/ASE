#DEP : addhead.pl (Certain bam files formed without header so added bthe SAM header manually)
#DEP : rm.multimappedreads.pl
###Test for 19 sam file
#created a errorlist of 19 files
# taking the the reads appeared more than twice and taking it in a order 

errordir='/proj/b2012046/rani/analysis/gsnap/bamsort/errortest'
cd $errordir
cat errorlist | xargs -I% echo "/proj/b2012046/rani/analysis/gsnap/bam/"%".fastq.bam " >errorbamfiles.list 
cat  errorbamfiles.list | xargs -I% echo "samtools view "%" | awk '{print \$1;}' | sort | uniq -c >"%.reads2nmapped >cmds.sh
sh cmds.sh

# taking out the number of uniquely mapped reads and non uniqly mapped reads
cat reads2nmapped.files | xargs -I% basename % | xargs -I% echo "awk '\$1 == 2 {print \$0;}'" ${errordir}/% " | wc -l >" %.n.uniq >cmds.sh
sh cmds.sh
cat reads2nmapped.files | xargs -I% basename % | xargs -I% echo "awk '\$1 != 2 {print \$0;}'" ${errordir}/% " | wc -l >" %.n.nonuniq >cmds2.sh
sh cmds2.sh


#Extract multimapped reads
ls -1 ${errordir}/*.reads2nmapped >reads2nmapped.files
cat reads2nmapped.files | xargs -I% basename % | xargs -I% echo "awk '\$1 != 2 {print \$0;}'" ${errordir}/% "| awk '{print \$2;}'>" %.nonuniq >cmds.sh
sh cmds.sh &

# Convert all bams to sams
bamdir='/proj/b2012046/rani/analysis/gsnap/bam'
sam=${errordir}
module load bioinfo-tools
module load samtools/0.1.18
cat  errorbamfiles.list | xargs -I% echo "samtools view -h -o /proj/b2012046/rani/analysis/gsnap/bamsort/errortest/"%".sam  "%".bam" >cmd.sh
(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2012046
#SBATCH -t 6:00:00
#SBATCH -J bamtosam
#SBATCH -p node -n 2
#SBATCH -e /proj/b2012046/rani/analysis/gsnap/bamsort/errortest/log/bamtosam.stderr
#SBATCH -o /proj/b2012046/rani/analysis/gsnap/bamsort/errortest/log/bamtosam.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=alva.rani@scilifelab.se
sh cmd.sh
EOF
)>${errordir}

#Rm multimapped reads
ls -1 *d.nonuniq | awk -F'.' -v OFS='.' '{print $1, $2;}' >errorfiles.prefix.list
cat errorfiles.prefix.list | xargs -I% echo "perl rm.multimappedreads.pl " %.sam %.fastq.bam.reads2nmapped.nonuniq ">"%.uniqmapped.sam >cmds.sh

# Convert all bams to sams
module load bioinfo-tools
module load samtools/0.1.18
ls  -1 *.uniqmapped.sam | xargs -I% echo  "samtools view -Sb -" >"%".bam >cmds.sh
(cat <<EOF
#!/bin/bash -l
#SBATCH -A b2012046
#SBATCH -t 6:00:00
#SBATCH -J samtobam
#SBATCH -p node -n 2
#SBATCH -e /proj/b2012046/rani/analysis/gsnap/bamsort/errortest/log/bamtosam.stderr
#SBATCH -o /proj/b2012046/rani/analysis/gsnap/bamsort/errortest/log/bamtosam.stdout
#SBATCH --mail-type=All
#SBATCH --mail-user=alva.rani@scilifelab.se
sh cmds.sh
EOF
)>${errordir}

# Additional steps
# Certain samples did'nt contain the header after the emoval of multi-mapped reads, for that sample I added header manually by
perl addheader.pl header 'The sample without header'


