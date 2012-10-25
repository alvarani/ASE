#median : #0.494
#mean:#0.495
sys = 'kalkyl'
sys = 'local'

if(sys == 'kalkyl'){
}
if(sys == 'local'){  
basedir ='/proj/b2012046/rani/analysis/degner'
  varcalls.datadir = file.path(basedir, 'varcall')
  basecount.datadir = file.path(basedir, 'basecount')  
  ase.datadir = file.path(basedir, 'ase')
  resdir = '/proj/b2012046/rani/analysis/degner/ase'
  annovar.datadir = file.path(basedir, 'annovar')
}

#mpileup calls
calls.regex = '.*hetvars.alleles$'

#het.vars
vars.regex = '.*hetvars$'

#basecounts
bc.regex = '.*basecount$'

#annovar
annovar.file = file.path(annovar.datadir, 'annovar_all_variants.RData')

#Res
bc.file = file.path(basecount.datadir, 'bc.RData')
ase.res.file = file.path(ase.datadir, 'ase.res.RData')
ase.sample.cond.diff.file = file.path(ase.datadir, 'induced.ase.res.RData')
induced.ase.sig.file = file.path(ase.datadir, 'induced.ase.sig.RData')
sig.ase.res.file = file.path(resdir, 'sig.ase.res.RData')

main <- function(){

  ###
  #Read data
  ###
  calls.files = list.files(varcalls.datadir, pattern = calls.regex, full.names = TRUE)
  vars.files = list.files(varcalls.datadir, pattern = vars.regex, full.names = TRUE)
  basecount.files = list.files(basecount.datadir, pattern = bc.regex, full.names = TRUE)

  samples = gsub('\\..*', '', basename(calls.files))
  samples = gsub('_12h', '', samples)
  basecount.list = list()
  length(basecount.list) = length(samples)
  names(basecount.list) = samples
  for(jsample in samples){
    basecount.file = basecount.files[grep(jsample, basecount.files)]
    calls.file = calls.files[grep(jsample, calls.files)]
    vars.file = vars.files[grep(jsample, vars.files)]
    basecount.list[[jsample]] = read.data(basecount.file, calls.file, vars.file)
  }

  #Dump
  save(basecount.list, file = bc.file)

  


  ###
  #Filter
  ###
  load(bc.file)
  bc.filt = basecount.list
    
  #Filter on DP
  min.dp = 10
  bc.filt = lapply(bc.filt, function(bc){bc[which(bc[, 'site.dp'] >= min.dp), ]})

  #total number of vars
  length(unique(unlist(lapply(bc.filt, '[[', 'site')))) #72377
  



  ###  
  #DP dist
  ###
  depths.pdf = file.path(resdir, sprintf('depths.mindp%i.pdf', min.dp))
  pdf(file = depths.pdf)
  plot.sample.dens(bc.filt, data.col = 'site.dp', main = sprintf('Read depth, minDP=%i', min.dp), xlab = 'dp', xlim = c(0, 500))
  dev.off()

  
  ###
  #Allele fracs dist (all vars)
  ###
  allelefracs.pdf = file.path(resdir, sprintf('allelefracs.mindp%i.pdf', min.dp))
  pdf(file = allelefracs.pdf)
  plot.sample.dens(bc.filt, data.col = 'alt.frac', main = sprintf('Alternate allele fraction, minDP=%i', min.dp), )
  abline(v=0.5)
  dev.off()


  #####
  #Binomial test
  #####
  #non-sim:
  #sig.vars = get.sig.mat(basecount.testres.list, alpha, pval.col = 'mod.padj')
  
  #Test
  ase.res = my.binom.test(bc.filt, alt.dp.col = 'alt.dp', site.dp.col = 'site.dp', p = 0.5)
  
  #Multiple testing correction
  ase.res = lapply(ase.res, function(sample.vars){padj = p.adjust(sample.vars[, 'pval'], method = 'BH'); sample.vars = cbind(sample.vars, padj); return(sample.vars);})

  #Order
  ase.res = lapply(ase.res, function(sample.vars){sample.vars = sample.vars[order(sample.vars[, 'padj']), ]; return(sample.vars);})

  #Dump
  save(ase.res, file = ase.res.file)

  
  ###
  #Sig hits stats
  ###
  load(ase.res.file)
  alpha = 0.05
  pval.col = 'padj'
  
  #Get sig hits
  sig.vars = lapply(ase.res, function(sample.vars, alpha, pval.col){sample.vars = sample.vars[which(sample.vars[, pval.col] <= alpha), ]; return(sample.vars);}, alpha = alpha, pval.col = pval.col)

 ## For all alternative variants 
  alt.vars = unlist(alt.vars)
  length(alt.vars) # 9299

  #Filter on alternative allele direction (significant variants)
  alt.vars = lapply(sig.vars, alt.filter)
  alt.sig.vars = unique(unlist(alt.vars))
  length(alt.sig.vars) # 8555
  
  #Dump
  save(sig.vars, alt.sig.vars, file = sig.ase.res.file)

  
  ###
  #Allele fracs dist (sig vars)
  ###
  allelefracs.pdf = file.path(resdir, sprintf('sig.allelefracs.mindp%i.pdf', min.dp))
  pdf(file = allelefracs.pdf)
  plot.sample.dens(sig.vars, data.col = 'alt.frac', main = sprintf('Alternate allele fraction, minDP=%i', min.dp), )
  abline(v=0.5)
  dev.off()
  
  #number of uniq sig.vars
  length(unique(unlist(lapply(sig.vars, '[[', 'site')))) #17107


    #number of uniq sig.vars, without filtering 
   length(unique(unlist(lapply(sig.vars, '[[', 'site'))))
  ### from here i havent run towards top:^

 
  
  ###
  #Gene based analysis. 
  ###
  #See: multiplesnps2gene.R.

} ## close for main function....

alt.filter <- function(jsample.ase, site.col = 'site', frac.col = 'alt.frac'){

  alt.higher.snvs = jsample.ase[which(jsample.ase[, frac.col] >0.5), site.col]
  
  return(alt.higher.snvs)
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

get.genomefrac <- function(hap, pos){
#hap.A was treated as the alternative allele when alt.frac was calculated
  
  a = merge(hap, pos, by = 'site')
  
  n.vars = nrow(a)
  A.genref.ind = which(a[, 'ref.A'] == a[, 'genome.ref'])
  B.genref.ind = which(a[, 'ref.B'] == a[, 'genome.ref'])  
  genome.alt.frac = vector(mode = 'numeric', length = n.vars)
  genome.alt.frac[A.genref.ind] = 1 - a[A.genref.ind, 'alt.frac']
  genome.alt.frac[B.genref.ind] = a[B.genref.ind, 'alt.frac']

  a = cbind(a, genome.alt.frac)
  return(a)
}

get.hapmatrix <- function(tab, id.col = 'feature_id'){

  mmseq = tab
  
  #get haplotype specific transcripts ids
  #(ALT: add ENST column with stripped haplo suffix, then tapply on ENST col)
  hapA.enst = unique(mmseq[grep('.*_A', mmseq[, id.col]), id.col])
  hapB.enst = unique(mmseq[grep('.*_B', mmseq[, id.col]), id.col])

  hapA.mmseq = merge(mmseq, as.matrix(hapA.enst), by.x = id.col, by.y = 1)
  hapB.mmseq = merge(mmseq, as.matrix(hapB.enst), by.x = id.col, by.y = 1)

  colnames(hapA.mmseq) = paste(colnames(hapA.mmseq), 'A', sep = '.')
  colnames(hapB.mmseq) = paste(colnames(hapB.mmseq), 'B', sep = '.')  

  #merge haps into one matrix
  nonhap.id.col = gsub('_A', '', hapA.mmseq[, paste(id.col, 'A', sep='.')])
  hapA.mmseq = cbind(hapA.mmseq, nonhap.id.col, stringsAsFactors = FALSE)
  nonhap.id.col = gsub('_B', '', hapB.mmseq[, paste(id.col, 'B', sep = '.')])
  hapB.mmseq = cbind(hapB.mmseq, nonhap.id.col, stringsAsFactors = FALSE)
  
  hap.mmseq = merge(hapA.mmseq, hapB.mmseq, by = 'nonhap.id.col')
  rownames(hap.mmseq) = hap.mmseq[, 'nonhap.id.col']
  
  #Rm redundant cols
  hap.mmseq = hap.mmseq[, setdiff(colnames(hap.mmseq), colnames(hap.mmseq)[grep('nonhap.id.col', colnames(hap.mmseq))])]
  hap.mmseq = hap.mmseq[, setdiff(colnames(hap.mmseq), colnames(hap.mmseq)[grep(id.col, colnames(hap.mmseq))])]
  
  return(hap.mmseq)
}


vars.list2mat <- function(vars, identical.cols = c('site', 'obs', 'ref')){
  samples = names(vars)
  n.samples = length(samples)
  all.sites = unique(unlist(lapply(vars, rownames)))
  n.sites = length(all.sites)
  sample.colnames = setdiff(colnames(vars[[1]]), identical.cols)
  samples.colnames = vector()
  for(jsample in samples){
    samples.colnames = c(samples.colnames, paste(jsample, sample.colnames, sep = '.'))
  }
  samples.colnames = c(identical.cols, samples.colnames)
  vars.mat = as.data.frame(matrix(NA, nrow = n.sites, ncol = length(samples.colnames), dimnames = list(all.sites, samples.colnames)))
  for(jsample in samples){
    print(jsample)
    var.sites = rownames(vars[[jsample]])    
    sample.cols = c(identical.cols, colnames(vars.mat)[grep(jsample, colnames(vars.mat))])
    sample.cols.noprefix = gsub(paste(jsample, '.', sep=''), '', sample.cols)
    vars.mat[var.sites, sample.cols] = vars[[jsample]][, sample.cols.noprefix]
  }  
  return(vars.mat)
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

plot.sample.dens <- function(data.list, data.col, main = 'allele frac', xlab = 'Allele fraction', xlim = NA){

  samples = names(data.list)
  density.list = list()
  for(jsample in 1:length(samples)){
    density.list[[jsample]] = density(data.list[[jsample]][, data.col])
  }  
  ylim = range(unlist(lapply(density.list, '[', 'y')))
  ylim[1] = floor(ylim[1])
  ylim[2] = ylim[2] + 0.1 * ylim[2]
  if(is.na(xlim)){
    xlim = range(unlist(lapply(density.list, '[', 'x')))
  }
  
  samples = names(data.list)
  for(jsample in 1:length(samples)){
    if(jsample == 1){
      plot(density.list[[jsample]], col = jsample, main = main, xlab = xlab, ylim = ylim, xlim = xlim)
    }
    else{
      lines(density.list[[jsample]], col = jsample, ylim = ylim, xlim = xlim)
    }
  }
  legend('topright', legend = samples, col = 1:length(samples), lty = 1, cex = 0.8)
}

read.data <- function(basecount.file, varcalls.file, hetvars.file){
  
  #Read basecounts
  basecount = read.table(basecount.file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
  colnames(basecount) = c('chr', 'pos', 'ref', 'alt', 'site.dp', 'ref.dp', 'alt.dp', 'alt.frac')
  site = paste(basecount[, 'chr'], basecount[, 'pos'], sep = '.')
  basecount = cbind(site, basecount, stringsAsFactors = FALSE)
  #Rm non-uniqe pos (ie, non-diallelic)
  site.n = table(site)
  diallelic.sites = names(site.n)[which(site.n == 1)]
  basecount = merge(basecount, as.matrix(diallelic.sites), by.x = 'site', by.y = 1)
  rownames(basecount) = basecount[, 'site']

  #Read varcalls (basecounts also contains an "alt" field, but if using I16, then it often includes undetermined alleles, "X")
  varcalls = read.table(varcalls.file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
  colnames(varcalls) = c('chr', 'pos', 'ref', 'alt')
  site = paste(varcalls[, 'chr'], varcalls[, 'pos'], sep = '.')
  varcalls = cbind(site, varcalls, stringsAsFactors = FALSE)
  #Rm non-uniqe pos (ie, non-diallelic)
  site.n = table(site)
  diallelic.sites = names(site.n)[which(site.n == 1)]
  varcalls = merge(varcalls, as.matrix(diallelic.sites), by.x = 'site', by.y = 1)
  rownames(varcalls) = varcalls[, 'site']  
  
  #Merge varcalls and basecounts
  change.col = c('alt')
  colnames(basecount)[grep(paste('^', change.col, '$', sep=''), colnames(basecount))] = paste('basecount', change.col, sep = '.')
  change.cols = c('alt')
  colnames(varcalls)[grep(paste('^', change.col, '$', sep=''), colnames(varcalls))] = paste('mpileup', change.col, sep = '.')  
  basecount = merge(varcalls[, c('site', 'mpileup.alt')], basecount, by = 'site')
  rownames(basecount) = basecount[, 'site']
  
  #Filter on heterozygosity
  het.vars = read.table(hetvars.file, stringsAsFactors = FALSE)
  colnames(het.vars) = c('chr', 'pos')
  site = paste(het.vars[, 'chr'], het.vars[, 'pos'], sep = '.')

  basecount = basecount[intersect(site, basecount[, 'site']), ]
  
  return(basecount)
}
