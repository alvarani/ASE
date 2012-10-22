

sys = 'local.work'
sys = 'local'

if(sys == 'local'){
  ase.datadir = '/proj/b2012046/rani/analysis/gsnap/ase'
  annovar.datadir = '/proj/b2012046/rani/analysis/gsnap/annovar'  
}

ase.res.file = file.path(ase.datadir, 'ase.res.RData') #See: snv.ase.R
annovar.file = file.path(annovar.datadir, 'annovar_all_variants.RData')
vars.annot.list.file = file.path(ase.datadir, 'mpileup.vars.annot.RData')
genes.pass.file = file.path(ase.datadir, 'multisnp.genes.pass.RData')

main <- function(){

  #get vars with ase pvals
  load(ase.res.file) #ase.res

  #total number of vars
  length(unique(unlist(lapply(ase.res, '[[', 'site')))) #97630

  #add gene annot    
  load(annovar.file) #varcalls
  chrpos = gsub(' ', '', apply(varcalls[, c('chr', 'start')], 1, paste,  collapse = '.'))
  
  chrpos=as.data.frame(chrpos)
  chrpos$chrpos2=paste('chr', chrpos$chrpos,sep="")
  chrpos=chrpos[,2, drop=FALSE]
  colnames(chrpos)="chrpos"


  varcalls = cbind(varcalls, chrpos, stringsAsFactors = FALSE)


  vars.list = lapply(ase.res, function(sample.vars, vars.annot){	
	vars = merge(sample.vars, vars.annot, by.x = 'site', by.y = 'chrpos'); return(vars)}, vars.annot = varcalls)


  #Check that same number of vars before and after annot
  a = unlist(lapply(ase.res, nrow))
  b = unlist(lapply(vars.list, nrow))
  print(a)
  print(b)
  length(unique(unlist(lapply(ase.res, rownames)))) #97630
  nrow(varcalls) #99104
  #Looks quite bad, but leave for now.
  
  #Dump
  save(vars.list, file = vars.annot.list.file)
  

  #Get genes with multiple significant ase:s (min 2 het.vars, see get.gene.asefrac)
  load(vars.annot.list.file)
  feat.col = 'gene annot'
  pval.col = 'padj'
  alpha = 0.05
  genes2asefrac.list = lapply(vars.list, get.gene.asefrac, alpha = alpha, feat.col = feat.col, pval.col = pval.col)  

  #Filter on frac of sig ase:s in gene
  frac = 1
  genes2asefrac.filt = lapply(genes2asefrac.list, gene2asefrac.filter, frac = frac)
  lapply(genes2asefrac.filt, nrow)
  length(unique(unlist(lapply(genes2asefrac.filt, rownames)))) #2407

  #Filter on number of significant alt alleles per gene
  frac.col = 'alt.frac'
  pval.col = 'padj'
  sample2genes.altallele = lapply(vars.list, filter.altallele, frac.col = frac.col, pval.col = pval.col)
  sample2genes.asefrac = lapply(genes2asefrac.filt, rownames)  
  samples = names(sample2genes.altallele)
  sample2genes = lapply(samples, function(jsample, sample2genes.asefrac, sample2genes.altallele){intersect(sample2genes.asefrac[[jsample]], sample2genes.altallele[[jsample]])}, sample2genes.asefrac = sample2genes.asefrac, sample2genes.altallele = sample2genes.altallele)
  names(sample2genes) = samples

  #Get common genes across samples
  gene2nsamples = table(unlist(sample2genes))
  gene2nsamples = gene2nsamples[order(gene2nsamples, decreasing = TRUE)] 
  length(gene2nsamples) #1508
  
  #Filter on n.samples
  n.samples = 2
  gene2nsamples.filt = gene2nsamples[which(gene2nsamples >= n.samples)]
  sample2genes.filt = lapply(sample2genes, function(jsample2genes, pass.genes){print(pass.genes); print(jsample2genes); return(intersect(jsample2genes, pass.genes))}, pass.genes = names(gene2nsamples.filt))
  lapply(sample2genes.filt, length)
  genes.pass = unique(unlist(sample2genes.filt))
  length(genes.pass) #156

  save(genes.pass, gene2nsamples, file = genes.pass.file)
  
  #FDR: 
  #See fdr.R  
  
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

filter.altallele <- function(vars, frac.thr = 0.5, gene.col = 'gene annot', frac.col = 'frac', pval.col = 'mod.padj', min.alt = 1, alpha = 0.05){
  
  sig.alt.vars = vars[which(vars[, pval.col] < alpha & vars[, frac.col] > frac.thr), ]
  gene2nvars = table(sig.alt.vars[, gene.col])
  genes.pass = names(gene2nvars)[which(gene2nvars >= min.alt)]

  #vars = merge(vars, as.matrix(genes.pass), by.x = gene.col, by.y = 1)
  
  return(genes.pass)
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
