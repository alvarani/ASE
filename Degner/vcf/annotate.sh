

#####
#Annotate variants
#####
vcfdir='/proj/b2012046/rani/analysis/degner/varcall'
annovardir='/proj/b2012046/rani/analysis/degner/annovar'
varsfile=${annovardir}'/all.vars.annovar'

#Format for annovar
cd $vcfdir
cat *.hetvars.alleles | sort -u | awk -F'\t' -v OFS='\t' '{print $1, $2, $2, $3, $4;}' >$varsfile
wc -l $varsfile #72429 
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
scriptdir='/proj/b2012046/rani/scripts/degner/annovar'
perl ${exedir}/annovar_annot_singlefileDe.pl $projid $annovardir $scriptdir 'all.vars.annovar'
find $scriptdir -name '*.sbatch' | grep 'ensg' | xargs -I% sbatch %
find $scriptdir -name '*.sbatch' | grep -v 'ensg' | xargs -I% sbatch %

### Submitted ....9-10-2012


#Merge annot files

exedir='/bubo/home/h24/alvaj/glob/code/ASE/vcf'
sh ${exedir}/varcall_merge_post_annovar_master_wf.0.degner.sh

#Convert annot file to RData structure
annovardir='/proj/b2012046/rani/analysis/degner/annovar'
exedir='/bubo/home/h24/alvaj/glob/code/ASE/vcf'
Rscript ${exedir}/annovar_flat2rdata.R ${annovardir}'/annovar_all_variants.txt' 'annovar_all_variants' &
