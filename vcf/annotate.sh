

#####
#Annotate variants
#####
vcfdir='/proj/b2012046/edsgard/ase/sim/data/varcalls'
annovardir='/proj/b2011075/analysis/sim/data/annovar'
varsfile=${annovardir}'/all.vars.annovar'

#Format for annovar
cd $vcfdir
cat *.hetvars.alleles | sort -u | awk -F'\t' -v OFS='\t' '{print $1, $2, $2, $3, $4;}' >$varsfile
wc -l $varsfile #
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
exedir='/bubo/home/h26/edsgard/glob/code/ase'
scriptdir='/proj/b20110FIX/analysis/FIX/scripts/annovar'
perl ${exedir}/annovar_annot_singlefile.pl $projid $annovardir $scriptdir 'all.vars.annovar'
find $scriptdir -name '*.sbatch' | grep 'ensg' | xargs -I% sbatch %
find $scriptdir -name '*.sbatch' | grep -v 'ensg' | xargs -I% sbatch %

#Merge annot files
#TBD: Need to change dirs within this script!!!!
sbatch ${scriptdir}/varcall_merge_post_annovar_master_wf.0.sh

#Convert annot file to RData structure
annovardir='/Users/edsgard/Documents/postdoc/ase/sim/annovar'
execdir='/Users/edsgard/Dropbox/postdoc/projects/ase/code/vcf'
Rscript ${execdir}/annovar_flat2rdata.R ${annovardir}'/annovar_master_all_subject_variants.txt' 'annovar_all_variants' &
