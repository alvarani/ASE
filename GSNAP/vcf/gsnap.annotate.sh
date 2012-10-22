

#####
#Annotate variants
#####
vcfdir='/proj/b2012046/rani/analysis/gsnap/varcalls'
annovardir='/proj/b2012046/rani/analysis/gsnap/annovar'
varsfile=${annovardir}'/all.vars.annovar'

#Format for annovar
cd $vcfdir
cat *.hetvars.alleles | sort -u | awk -F'\t' -v OFS='\t' '{print $1, $2, $2, $3, $4;}' >$varsfile
wc -l $varsfile #878
cat $varsfile | awk '{print $1;}' | sort -u >chroms

#Split into chromwise files
chroms=($(cat ${vcfdir}/chroms))
for chr in ${chroms[@]}
do
    echo ${chr}
    awk -F'\t' -v OFS='\t' -v chrom=${chr} '$1 == chrom' ${varsfile} >${varsfile}.${chr}
done

#Annotate

projid='b2012046'
exedir='/bubo/home/h24/alvaj/glob/code/ASE/vcf'
scriptdir='/proj/b2012046/rani/scripts/gsnap/annovar'
perl ${exedir}/annovar_annot_singlefile.pl $projid $annovardir $scriptdir 'all.vars.annovar'
find $scriptdir -name '*.sbatch' | grep 'ensg' | xargs -I% sbatch % ## need to submit this first and see the results.
find $scriptdir -name '*.sbatch' | grep -v 'ensg' | xargs -I% sbatch %
### Submitted .... 9-10-2012


#Merge annot files
#TBD: Need to change dirs within this script!!!!
exedir='/bubo/home/h24/alvaj/glob/code/ASE/vcf'
sh ${exedir}/varcall_merge_post_annovar_master_wf.0.sh

#Convert annot file to RData structure
annovardir='/proj/b2012046/rani/analysis/gsnap/annovar'
exedir='/bubo/home/h24/alvaj/glob/code/ASE/vcf'
Rscript ${exedir}/annovar_flat2rdata.R ${annovardir}'/annovar_all_variants.txt' 'annovar_all_variants' &

#my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT");
