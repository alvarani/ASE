
#USAGE: Rscript bedtoolshist2pdf.R histfilelist.file

#argv = commandArgs(trailingOnly=TRUE)
#histfilelist.file = argv

#annot = 'ccds'
histfilelist.file = '/proj/b2012046/rani/analysis/degner/ccds.histfiles'
resdir = dirname(histfilelist.file)
tab.outfile = file.path(resdir, 'ccds.cumsum.tab')
hist.pdf.file = file.path(resdir, 'ccds.hist.pdf')

#annot = 'genome'
#histfilelist.file = '../res/hs_pipe/run1/coverage/genome.histfiles'
#resdir = dirname(histfilelist.file)
#tab.outfile = file.path(resdir, 'genome.cumsum.tab')
#hist.pdf.file = file.path(resdir, 'genome.hist.pdf')


#read data
histfile2id = read.table(histfilelist.file, stringsAsFactors=FALSE, row.names=1)
colnames(histfile2id) = 'filename'
#hist.files = c('../res/hs_pipe/test/bedtools.genome.all.hist', '../res/hs_pipe/test/bedtools.ensembl.all.hist', '../res/hs_pipe/test/bedtools.ccds.all.hist')

#set linetype (dashed for unstim)
samples = rownames(histfile2id)
n.samples = length(samples)
ltypes = rep(1, n.samples)
names(ltypes) = samples
ltypes[grep('unstim', names(ltypes))] = 2

#set color
library('RColorBrewer')
individual = gsub('^(\\d+)_.*', '\\1', samples, perl = TRUE)
ind2color = brewer.pal(8, 'Dark2')
names(ind2color) = unique(individual)
colors = ind2color[individual]

histfile2id = cbind(histfile2id, ltypes, colors, stringsAsFactors=FALSE)

main <- function(histfile2id, tab.outfile, hist.pdf.file){
  hist.tabs = read.files(histfile2id)
  dp2frac = get.dp2frac(hist.tabs)
  dp2bases = get.dp2bases(hist.tabs)
  plot.dp2frac(dp2frac, ltypes = histfile2id[names(dp2frac), 'ltypes'], colors = histfile2id[names(dp2frac), 'colors'])
  plot.dp2bases(dp2bases)


  #Plot
  if(annot == 'ccds'){
    ylab = 'bases covered'
    pdf(file = hist.pdf.file, width = 10, height = 7)
    par(mfrow = c(1, 2))
    unstim.ind = grep('unstim', names(dp2bases))
    dp2bases.unstim = dp2bases[unstim.ind]
    plot.dp2y(dp2bases.unstim, ltypes = histfile2id[names(dp2bases.unstim), 'ltypes'], colors = histfile2id[names(dp2bases.unstim), 'colors'], ylab = ylab, min.plot.dp=2)
    
    lps.ind = grep('LPS', names(dp2bases))
    dp2bases.lps = dp2bases[lps.ind]
    plot.dp2y(dp2bases.lps, ltypes = histfile2id[names(dp2bases.lps), 'ltypes'], colors = histfile2id[names(dp2bases.lps), 'colors'], ylab = ylab, min.plot.dp = 2)
    dev.off()
  } else{
    #Plot
    pdf(file = hist.pdf.file, width = 10, height = 7)
    par(mfrow = c(1, 2))
    unstim.ind = grep('unstim', names(dp2frac))
    dp2frac.unstim = dp2frac[unstim.ind]
    plot.dp2frac(dp2frac.unstim, ltypes = histfile2id[names(dp2frac.unstim), 'ltypes'], colors = histfile2id[names(dp2frac.unstim), 'colors'])
    
    lps.ind = grep('LPS', names(dp2frac))
    dp2frac.lps = dp2frac[lps.ind]
    plot.dp2frac(dp2frac.lps, ltypes = histfile2id[names(dp2frac.lps), 'ltypes'], colors = histfile2id[names(dp2frac.lps), 'colors'])
    dev.off()
  }
  
  #Print table
  max.dp = 101
  dp2frac.mat = dp2frac[[1]][1:max.dp]
  for(jfile in 2:length(dp2frac)){
    dp2frac.mat = cbind(dp2frac.mat, dp2frac[[jfile]][1:max.dp])
  }
  colnames(dp2frac.mat) = names(dp2frac)
  
  rownames(dp2frac.mat) = 0:(max.dp - 1)
  write.table(dp2frac.mat, sep = '\t', quote = FALSE, col.names = NA, file = tab.outfile)

  #avg cov
  avg.cov = unlist(lapply(hist.tabs, function(hist.tab){sum(as.numeric(hist.tab[, 'dp'] * hist.tab[, 'bases.cov'])) / hist.tab[1, 'feat.len']}))  
  write.table(round(avg.cov), sep='\t', quote = F, col.names = F, file = file.path(resdir, 'ccds.avgcov.tab'))

  #y-axis: n.genes, x-axis: avg.cov
}

read.files <- function(histfile2id){

  samples = rownames(histfile2id)
  n.hist = length(samples)
  
  hist.tabs = list()
  length(hist.tabs) = n.hist
  names(hist.tabs) = samples
  for(jsample in samples){
    hist.tab = read.table(histfile2id[jsample, 'filename'], sep='\t', stringsAsFactors=FALSE)
    colnames(hist.tab) = c('feat', 'dp', 'bases.cov', 'feat.len', 'frac.cov')
    hist.tabs[[jsample]] = hist.tab
  }

  return(hist.tabs)
}

get.dp2frac <- function(hist.tabs){
#fraction of bases covered against dp

  n.hist = length(hist.tabs)
  dp2frac = list()
  length(dp2frac) = n.hist
  hist.files = names(hist.tabs)
  names(dp2frac) = hist.files
  for(hist.file in hist.files){
    hist.tab = hist.tabs[[hist.file]]
    dp2frac[[hist.file]] = rev(cumsum(hist.tab[rev(1:nrow(hist.tab)), 'frac.cov']))
  }

  return(dp2frac)
}

get.dp2bases <- function(hist.tabs){
#number of bases covered against dp
  n.hist = length(hist.tabs)
  dp2bases = list()
  length(dp2bases) = n.hist
  hist.files = names(hist.tabs)
  names(dp2bases) = hist.files
  for(hist.file in hist.files){
    hist.tab = hist.tabs[[hist.file]]
    dp2bases[[hist.file]] = rev(cumsum(hist.tab[rev(1:nrow(hist.tab)), 'bases.cov']))
  }

  return(dp2bases)
}

plot.dp2y <- function(dp2frac, max.plot.dp = 100, ylab = 'fraction covered', ltypes, colors, min.plot.dp = 1){

  n.hist = length(dp2frac)
  
  ylim = range(unlist(dp2frac))
  ylim = c(0, 5e8)
  jhist = 1
  plot(dp2frac[[jhist]][min.plot.dp:max.plot.dp], xlab = 'depth', ylab = ylab, col=colors[jhist], type = 'l', ylim = ylim, lty = ltypes[jhist])
  for(jhist in 2:n.hist){
    lines(dp2frac[[jhist]][min.plot.dp:max.plot.dp], col=colors[jhist], lty = ltypes[jhist])
  }
  legend('topright', legend = names(dp2frac), lty = ltypes, col = colors)
}

#execute
#main(histfile2id)
