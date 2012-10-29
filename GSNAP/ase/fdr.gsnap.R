

sys = 'kalk'
sys = 'local'

if(sys == 'kalk'){
  
                  }
if(sys == 'local'){
                   syntdir = '/proj/b2012046/edsgard/ase/sim/data/synt/ase'
                   ase.method.resdir = '/proj/b2012046/rani/analysis/gsnap/ase'

                 }
#Variant-based analysis
synt.sig.ase.res.file = file.path(syntdir, 'ase.noninduced.vars.RData')
sig.ase.res.file = file.path(ase.method.resdir, 'sig.ase.res.RData')
#Gene-based analysis
synt.genes.pass.file = file.path(syntdir, 'multisnp.genes.pass.RData')
method.genes.pass.file = file.path(ase.method.resdir, 'multisnp.genes.pass.RData')

#Cond-dep ASE
induced.ase.sig.file = file.path(syntdir, 'induced.ase.sig.RData')


main <- function(){

  
  ######
  #1. Condition-INDEPENDENT ASE
  ######

  ###
  #1.1 Variants
  ###

  #Load "TRUTH"
  load(synt.sig.ase.res.file) #sig.vars, alt.sig.vars
  true.sig.vars = unique(unlist(lapply(sig.vars, '[[', 'chrpos')))
  true.alt.sig.vars = alt.sig.vars

  #Load significant variants after application of ASE method
  load(sig.ase.res.file) #sig.vars, alt.sig.vars
  sig.vars = unique(unlist(lapply(sig.vars, rownames)))
  sig.vars = gsub('^chr', '', sig.vars)
  alt.sig.vars = gsub('^chr', '', alt.sig.vars)

  #Specificity
  #get fdr, alt allele direction filtered: NO
  fdr = get.fdr(sig.vars, true.sig.vars)
  print(fdr) #Specificity:7.4%

  #get fdr, alt allele direction filtered: YES
  fdr = get.fdr(alt.sig.vars, true.alt.sig.vars)
  print(fdr) #Specificity: 6.3%
  
  #SENSITIVITY
  #get sensitivity
  sens = get.sens(alt.sig.vars, true.alt.sig.vars)
  print(sens) #1.24 %
  # For signifivant variants
  sens = get.sens(sig.vars, true.alt.sig.vars)
  print(sens) #2.07%
  
  ###
  #1.2 Genes
  ###
  
  #Load "TRUTH"
  load(synt.genes.pass.file) #genes.pass, gene2nsamples
  true.genes.pass = genes.pass
  true.gene2nsamples = gene2nsamples

  #Load significant variants after application of ASE method
  load(method.genes.pass.file) #genes.pass, gene2nsamples
  
  #get fdr, n.samples = 2
  fdr = get.fdr(genes.pass, true.genes.pass)
  print(fdr) #46%, 72, 156

  #n.samples = 1
  true.genes.pass = names(true.gene2nsamples)
  genes.pass = names(gene2nsamples)
  fdr = get.fdr(genes.pass, true.genes.pass)
  print(fdr) #19%, 300, 1508

  #SENSITIVITY
  #get sensitivity
   #n.samples = 2
  sens = get.sens(genes.pass, true.genes.pass)
  print(sens) #1%  156.0000000 1352.0000000
  # For signifivant variants
  #n.samples = 1
  true.genes.pass = names(true.gene2nsamples)
  genes.pass = names(gene2nsamples)
  sens = get.sens(genes.pass, true.genes.pass)
  print(sens) #9.6%
    
#********# my results#**************

  ##########
  #2. Condition-dependent ASE
  ##########
  load(induced.ase.sig.file) #sig.vars, genes.pass, sig.vars.nsamples1, genes.pass.nsamples1
  true.sig.vars = sig.vars
  true.genes.pass = genes.pass
  true.sig.vars.nsamples1 = sig.vars.nsamples1
  true.genes.pass.nsamples1 = genes.pass.nsamples1
  induced.ase.sig.file = file.path(ase.method.resdir, 'induced.ase.sig.RData')
  load(induced.ase.sig.file) #sig.vars, genes.pass

  ##
  #2.1 Variants
  ##
  #1 individual
  fdr = get.fdr(rownames(sig.vars.nsamples1), rownames(true.sig.vars.nsamples1))
  print(fdr) #36%, 2137, 5865
  
  #>=2 individuals
  fdr = get.fdr(rownames(sig.vars), rownames(true.sig.vars))
  print(fdr) #44%, 435, 998. OLD: 56%, 557, 998

  ##
  #2.2 Genes
  ##
  #1 individual
  fdr = get.fdr(genes.pass.nsamples1, true.genes.pass.nsamples1)
  print(fdr) #7.5%, 69, 917
  
  #>=2 individuals
  fdr = get.fdr(genes.pass, true.genes.pass)
  print(fdr) #17%, 12, 69
  
}

alt.filter <- function(jsample.ase){
  b.ref.snvs = jsample.ase[which(jsample.ase[, 'ref.1kg'] == jsample.ase[, 'hapB.allele']), 'chrpos']
  a.ref.snvs = jsample.ase[setdiff(rownames(jsample.ase), b.ref.snvs), 'chrpos']
  
  alt.higher.snvs = union(b.ref.snvs[which(jsample.ase[b.ref.snvs, 'hapB.frac'] < 0.5)], a.ref.snvs[which(jsample.ase[a.ref.snvs, 'hapB.frac'] > 0.5)])
  return(alt.higher.snvs)
}

get.sig.mat <- function(all.vars, alpha, pval.col = 'mod.padj'){
  sig.vars = lapply(all.vars, function(jsample.vars, alpha){sig.vars = jsample.vars[which(jsample.vars[, pval.col] <= alpha), ]; return(rownames(sig.vars));}, alpha = alpha)
  
  samples = names(sig.vars)
  all.sig.vars = unique(unlist(sig.vars))
  sig.vars.mat = matrix(0, nrow=length(all.sig.vars), ncol = length(sig.vars), dimnames = list(all.sig.vars, samples))
  for(jsample in samples){
    sig.vars.mat[sig.vars[[jsample]], jsample] = 1
  }
  return(sig.vars.mat)
}

get.ase.sample.cond.diff <- function(sim.list, site.dp.col = 'site.dp', alt.dp.col = 'allele.dp'){
  
  samples = unique(gsub('_.*', '', names(sim.list)))

  ase.sample.cond.diff.list = list()
  length(ase.sample.cond.diff.list) = length(samples)
  names(ase.sample.cond.diff.list) = samples
  for(jsample in samples){
    print(jsample)
    untreat.vars = sim.list[[paste(jsample, 'unstim', sep = '_')]]
    treat.vars = sim.list[[paste(jsample, 'LPS', sep = '_')]]

    shared.vars = merge(treat.vars, untreat.vars, by = 'row.names') #ADDs Xs and Ys to colnames
    rownames(shared.vars) = shared.vars[, 'Row.names']
    shared.vars = shared.vars[, setdiff(colnames(shared.vars), 'Row.names')]
  
    #Fisher test
    counts = shared.vars[, c(paste(site.dp.col, 'x', sep = '.'), paste(alt.dp.col, 'x', sep = '.'), paste(site.dp.col, 'y', sep='.'), paste(alt.dp.col, 'y', sep = '.'))]
    ref.dp.x = counts[, paste(site.dp.col, 'x', sep = '.')] - counts[, paste(alt.dp.col, 'x', sep = '.')]
    ref.dp.y = counts[, paste(site.dp.col, 'y', sep = '.')] - counts[, paste(alt.dp.col, 'y', sep = '.')]
    counts = cbind(counts, ref.dp.x, ref.dp.y)
    ase.sample.cond.diff = t(apply(counts, 1, function(var.counts){mat = matrix(var.counts[c(paste(alt.dp.col, 'x', sep = '.'), 'ref.dp.x', paste(alt.dp.col, 'y', sep = '.'), 'ref.dp.y')], byrow=TRUE, nrow=2); res = fisher.test(mat); return(c(res$p.value, res$estimate, res$conf.int[1:2]));}))
    colnames(ase.sample.cond.diff) = c('pval', 'or', 'or.ci.lower', 'or.ci.upper')

    #multiple testing correction
    pval.bh = p.adjust(ase.sample.cond.diff[, 'pval'], method = 'BH')
    ase.sample.cond.diff = cbind(ase.sample.cond.diff, pval.bh)
    ase.sample.cond.diff = ase.sample.cond.diff[order(ase.sample.cond.diff[, 'pval']), ]

    #store
    ase.sample.cond.diff.list[[jsample]] = ase.sample.cond.diff
  }
  return(ase.sample.cond.diff.list)
}

get.fdr <- function(feat.pass, true.feat.pass){
  tp = intersect(feat.pass, true.feat.pass)
  fp = setdiff(feat.pass, tp)
  n.fp = length(fp)
  n.p = length(feat.pass)
  fdr = n.fp / n.p

  fdr = c(fdr, n.fp, n.p)
  return(fdr)
}

get.sens <- function(feat.pass, true.feat.pass){
  tp = intersect(feat.pass, true.feat.pass)
  fn = setdiff(true.feat.pass, feat.pass)
  n.tp = length(tp)
  n.fn = length(fn)
  sens = n.tp / (n.tp + n.fn)

  sens = c(sens, n.tp, n.fn)
  return(sens)
}

filter.altallele <- function(vars, frac.thr = 0.5, gene.col = 'gene annot', frac.col = 'frac', pval.col = 'mod.padj', min.alt = 1, alpha = 0.05){
  
  sig.alt.vars = vars[which(vars[, pval.col] < alpha & vars[, frac.col] > frac.thr), ]
  gene2nvars = table(sig.alt.vars[, gene.col])
  genes.pass = names(gene2nvars)[which(gene2nvars >= min.alt)]

  #vars = merge(vars, as.matrix(genes.pass), by.x = gene.col, by.y = 1)
  
  return(genes.pass)
}

gene2asefrac.filter <- function(genes2asefrac, frac){

  genes2asefrac.filt = genes2asefrac[which(genes2asefrac[, 'frac.sig'] >= frac), ]

  return(genes2asefrac.filt)
}

get.sigvars <- function(vars, alpha, pval.col = 'padj'){

  pval.cols = colnames(vars)[grep(pval.col, colnames(vars))]
  if(length(pval.cols) > 1){
    sig.vars = names(which(apply(vars[, pval.cols], 1, function(pvals, alpha){any(pvals < alpha, na.rm = TRUE)}, alpha = alpha)))
    sig.vars = vars[sig.vars, ]
  }
  else{
    sig.vars = vars[which(vars[, pval.cols] < alpha), ]
  }

  return(sig.vars)
}

get.gene.asefrac <- function(vars, alpha = alpha, feat.col = 'hugo.enst', pval.col = 'padj'){
  
  #Count number of (tested) het.vars per gene
  gene2nvars = table(vars[, feat.col])
  gene2nvars = as.matrix(gene2nvars[order(gene2nvars, decreasing = TRUE)])
  colnames(gene2nvars) = 'n.het.vars'
  
  #Count number of sig.vars per gene
  sig.vars = get.sigvars(vars, alpha, pval.col = pval.col)
  gene2nsigvars = table(sig.vars[, feat.col])
  gene2nsigvars = as.matrix(gene2nsigvars[order(gene2nsigvars, decreasing = TRUE)])
  colnames(gene2nsigvars) = 'n.sig.vars'

  #Get frac of vars
  gene2asevars = merge(gene2nvars, gene2nsigvars, by = 'row.names')
  rownames(gene2asevars) = gene2asevars[, 'Row.names']
  gene2asevars = gene2asevars[, setdiff(colnames(gene2asevars), 'Row.names')]
  frac.sig = gene2asevars[, 'n.sig.vars'] / gene2asevars[, 'n.het.vars']
  gene2asevars = cbind(gene2asevars, frac.sig)

  #filter genes with at least 2 vars
  min.vars = 2
  gene2asevars.min2 = gene2asevars[which(gene2asevars[, 'n.het.vars'] >= min.vars), ]

  #order
  gene2asevars.min2 = gene2asevars.min2[order(gene2asevars.min2[, 'frac.sig'], gene2asevars.min2[, 'n.het.vars'], decreasing = TRUE), ]

  return(gene2asevars.min2)
}

save.annot <- function(annovar.tab.file, annovar.file){
  vars.annot = read.table(annovar.tab.file, stringsAsFactors = FALSE)
  vars.annot = vars.annot[, c(2,3,4)]
  colnames(vars.annot) = c('gene annot', 'chr', 'pos')
  chrpos = gsub(' ', '', apply(vars.annot[, c('chr', 'pos')], 1, paste, collapse = '.'))
  vars.annot = cbind(vars.annot, chrpos, stringsAsFactors = FALSE)
  rownames(vars.annot) = vars.annot[, 'chrpos']
  save(vars.annot, file = annovar.file)
}

get.allelefracs <- function(snp2count){
  
  snp2count.red = unique(snp2count[, c('chrpos', 'chr', 'pos', 'ref.1kg', 'alt.1kg', 'maf', 'hapA.allele', 'hapB.allele')])

  #Error check
  n.snps = nrow(snp2count.red)
  n.uniq.snps = length(unique(snp2count[, 'chrpos']))
  if(n.snps != n.uniq.snps){
    warning('There are several genotypes at single positions!')
      snp2count.red.list = split(snp2count.red, snp2count.red[, 'chrpos'])
    nonuniq.snps.ind = which(unlist(lapply(snp2count.red.list, nrow)) != 1)
    uniq.snps.ind = setdiff(1:length(snp2count.red.list), nonuniq.snps.ind)
    uniq.snps = names(snp2count.red.list)[uniq.snps.ind]
    snp2count.red.uniq = merge(snp2count.red, as.matrix(uniq.snps), by.x = 'chrpos', by.y = 1)
    rownames(snp2count.red.uniq) = snp2count.red.uniq[, 'chrpos']
  }

  #Sum counts for all ENST where SNV is found
  hapA.count = tapply(snp2count[, 'hapA.count'], snp2count[, 'chrpos'], sum)
  hapB.count = tapply(snp2count[, 'hapB.count'], snp2count[, 'chrpos'], sum)
  hapA.count = hapA.count[snp2count.red[, 'chrpos']]
  hapB.count = hapB.count[snp2count.red[, 'chrpos']]
  snp2count.red.uniq = cbind(snp2count.red, hapA.count, hapB.count, stringsAsFactors = FALSE)

  #Add site.dp and alt.frac cols
  site.dp = apply(snp2count.red.uniq[, c('hapA.count', 'hapB.count')], 1, sum)
  hapB.frac = snp2count.red.uniq[, 'hapB.count'] / site.dp
  snp2count.red.uniq = cbind(snp2count.red.uniq, site.dp, hapB.frac, stringsAsFactors = FALSE)

  return(snp2count.red.uniq)
}

my.binom.test <- function(vars, alt.dp.col = 'allele.dp', site.dp.col = 'site.dp', p = 0.5){
  samples = names(vars)
  for(jsample in samples){
    print(jsample)
    jsample.vars = vars[[jsample]]
    n.vars = nrow(jsample.vars)
    b.res = matrix(nrow = n.vars, ncol = 3)
    colnames(b.res) = c('pval', 'ci.lower', 'ci.upper')
    for(jvar in 1:n.vars){
      binom.res = binom.test(jsample.vars[jvar, alt.dp.col], jsample.vars[jvar, site.dp.col], p = p, alternative = 'two.sided')
      b.res[jvar, ] = c(binom.res$p.value, binom.res$conf.int[1:2])
    }
    vars[[jsample]] = cbind(jsample.vars, b.res)
  }
  return(vars)
}

read.data <- function(exprdir, snpdir, number2sample.map.file){
  #snp2enst. key: snp, values: ensts
  #enst2counts. key: enst, values: expr
  #merge: snp2enst, enst2count => snp2count
  
  #read enstcounts
  enstcounts.files = list.files(exprdir, 'sim.enst.counts.tab\\..*', full.names = TRUE)
  enst2counts.list = lapply(enstcounts.files, function(jfile){enst2counts = read.table(jfile, sep = '\t', stringsAsFactors = FALSE); colnames(enst2counts) = c('enst', 'tot', 'hapA', 'hapB'); return(enst2counts);})
  names(enst2counts.list) = basename(enstcounts.files)

  #read snps
  snp.files = list.files(snpdir, 'synt.snps.tab\\..*', full.names = TRUE)
  snps.list = lapply(snp.files, function(jfile){snps = read.table(jfile, stringsAsFactors = FALSE, sep = '\t'); colnames(snps) = c('enst', 'chr', 'pos', 'ref.1kg', 'alt.1kg', 'maf', 'hapA', 'hapB'); return(snps);})
  names(snps.list) = basename(snp.files)

  #map names
  number2sample.map = read.table(number2sample.map.file, sep = '\t', stringsAsFactors = FALSE)
  colnames(number2sample.map) = c('snpfile', 'enstfile')
  number2sample.map[, 'snpfile'] = gsub('customref', 'synt.snps.tab', number2sample.map[, 'snpfile'])
  number2sample.map[, 'snpfile'] = gsub('.fa.haplo.enst', '', number2sample.map[, 'snpfile'])
  sample = gsub('sim.enst.counts.tab.', '', number2sample.map[, 'enstfile'])
  number2sample.map = cbind(number2sample.map, sample)
  rownames(number2sample.map) = number2sample.map[, 'snpfile']
  
  names(snps.list) = number2sample.map[names(snps.list), 'sample']
  rownames(number2sample.map) = number2sample.map[, 'enstfile']
  names(enst2counts.list) = number2sample.map[names(enst2counts.list), 'sample']
  
  #merge
  samples = names(snps.list)
  snp2count.list = lapply(samples, function(jsample, snps.list, enst2counts.list){snp2count = merge(snps.list[[jsample]], enst2counts.list[[jsample]], by = 'enst'); colnames(snp2count) = gsub('x', 'allele', colnames(snp2count)); colnames(snp2count) = gsub('y', 'count', colnames(snp2count)); return(snp2count);}, snps.list = snps.list, enst2counts.list = enst2counts.list)
  names(snp2count.list) = samples

  #add site col
  snp2count.list = lapply(snp2count.list, function(snps){chrpos = gsub(' ', '', apply(snps[, c('chr', 'pos')], 1, paste, collapse = '.')); snps = cbind(snps, chrpos, stringsAsFactors = FALSE); return(snps);})
    
  return(snp2count.list)

  #check
  #lapply(snps.list, nrow)
  #lapply(snp2count.list, nrow)
  #site.list = lapply(snp2count.list, '[[', 'chrpos')
  #lapply(site.list, function(x){length(unique(x))})
  #Since many snps maps to several ENST...
}

old.check <- function(){

  
  #Accomodate the uniq/non-uniq issue
  load(snp2count.uniq.file) #snp2count.uniq.list, nonuniq.sites.list
  load(annovar.file) #vars.annot
  
  #map to gene annot
  samples = names(nonuniq.sites.list)
  nonuniq.genes.list = list()
  length(nonuniq.genes.list) = length(samples)
  names(nonuniq.genes.list) = samples
  for(jsample in samples){
    print(jsample)
    sample.vars = nonuniq.sites.list[[jsample]]
    nonuniq.sites = intersect(sample.vars, rownames(vars.annot))
    nonuniq.genes.list[[jsample]] = vars.annot[nonuniq.sites, 'gene annot'];
  }
  nonuniq.genes.list = lapply(nonuniq.sites.list, function(sample.vars, vars.annot){genes = vars.annot[intersect(sample.vars, rownames(vars.annot)), 'gene annot']; return(genes)}, vars.annot = vars.annot)

  #dump
  save(nonuniq.genes.list, file = nonuniq.genes.file)
  
  #rm these nonuniq.genes from genes.pass
  load(nonuniq.genes.file) #nonuniq.genes.list
  nonuniq.genes = unique(unlist(nonuniq.genes.list))
  length(nonuniq.genes) #15455
  genes.pass.uniq = setdiff(genes.pass, nonuniq.genes) #373 -> 31
  
  fdr = get.fdr(genes.pass.uniq, true.genes.pass)
  print(fdr) #52%
  #TBD: but which true.genes that passed is dep on the multiple testing correction. In the MTC also nonuniq vars were present...

  
  ###
  #Condition-Dependent ASE
  ###
  #rm these nonuniq.genes from genes.pass
  load(nonuniq.genes.file) #nonuniq.genes.list
  nonuniq.genes = unique(unlist(nonuniq.genes.list))
  length(nonuniq.genes) #15455
  genes.pass.uniq = setdiff(genes.pass, nonuniq.genes) #69 -> 2
  
  fdr = get.fdr(genes.pass.uniq, true.genes.pass)
  print(fdr) #100%. But only 2 genes!!
  #TBD: but which true.genes that passed is dep on the multiple testing correction. In the MTC also nonuniq vars were present...

  
  ##
  #Annot
  ##
  #just used refgene annot to avoid waiting for annovar annots to finish
  if(0){
    save.annot(annovar.tab.file, annovar.file)
  }

  
  #OLD: Need to dump the nonuniq.snps to check if among the identified hits..
  sites.list = lapply(snp2count.list, '[[', 'chrpos')
  uniq.sites.list = lapply(snp2count.uniq.list, '[[', 'chrpos')
  nonuniq.sites.list = lapply(samples, function(jsample, sites.list, uniq.sites.list){setdiff(sites.list[[jsample]], uniq.sites.list[[jsample]])}, sites.list = sites.list, uniq.sites.list = uniq.sites.list)
  names(nonuniq.sites.list) = samples

  ###
  #Check disconcordant vars...
  ###
  snp2count.red = unique(snp2count[, c('chrpos', 'chr', 'pos', 'ref.1kg', 'alt.1kg', 'maf', 'hapA.allele', 'hapB.allele')])
  hapA.count = tapply(snp2count[, 'hapA.count'], snp2count[, 'chrpos'], sum)
  length(hapA.count) #91363
  nrow(snp2count.red) #109882
  snp2count.red.list = split(snp2count.red, snp2count.red[, 'chrpos'])
  nonuniq.snps.ind = which(unlist(lapply(snp2count.red.list, nrow)) != 1)
  snp2count.red.list[head(nonuniq.snps.ind, 2)]
  length(nonuniq.snps.ind) #15824
  length(uniq.snps.ind) #75539

  
  #TBD: NB: Backtracing to the root of this...
  #grep '100336361' /Users/danieledsgard/Documents/postdoc/ase/sim/snps/synt.snps.tab.*
  #Clearly the problem is present in these files... That basically means that one single individual will not have a well-defined var such single positions:(
  #See makesyndata.R (reading synt.snps.tab, which contains several SNPs per ENST)

  #TBD: consider these as noise?
  #ALT A) Exclude them
  #ALT B) Don't distinguish AB from BA. Would then "just" need to swap alleles. However, the synthetic reads will be dependent on the actual haplotype... Could "just" exclude reads from conflicting enst from the fastq file.
  #ALT C) Keep haplo-info.

}
