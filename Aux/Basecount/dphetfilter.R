
#USAGE: Rscript dp.file gt.file 10

#get args
argv = commandArgs(trailingOnly=TRUE)
dp.file = argv[1]
gt.file = argv[2]
dp.min = as.integer(argv[3])

#Test-suite
if(0){
  datadir = '~/Dropbox/postdoc/projects/ase/data/codetest'
  dp.file = file.path(datadir, 'head.dp')
  gt.file = file.path(datadir, 'head.gt')
  dp.min = 10
}

main <- function(dp.file, gt.file, dp.min){
  
  #read
  dp.mat = read.dp(dp.file)
  gt.mat = read.gt(gt.file)

  #sort
  dp.mat = dp.mat[rownames(gt.mat), ]

  #filter
  rseq.pos.pass = dp.het.filter(dp.mat, gt.mat, dp.min)

  #dump
  samples = names(rseq.pos.pass)
  for(jsample in samples){
    hetvars.file = paste(jsample, dp.file, 'het', 'pos', sep='.')
    write.table(rseq.pos.pass[[jsample]], quote = F, sep = '\t', col.names=F, row.names = F, file = hetvars.file)
  }
}

dp.het.filter <- function(dp.mat, gt.mat, dp.min){
  
  dp.pass = lapply(dp.mat, function(ind.vals, minval){which(ind.vals >= minval)}, minval = dp.min)
  het.pass = lapply(gt.mat, function(ind.vals, filt.val){which(ind.vals == filt.val)}, filt.val = 1)
  samples = names(dp.pass)
  sample2pass = lapply(samples, function(sample, dp.pass, het.pass){intersect(dp.pass[[sample]], het.pass[[sample]]);}, dp.pass = dp.pass, het.pass = het.pass)
  rseq.pos.pass = lapply(sample2pass, function(ind.pass, pos){pos[ind.pass]}, pos = rownames(dp.mat))
  names(rseq.pos.pass) = samples
  
  return(rseq.pos.pass)
}

read.gt <- function(gt.file){

  gt.mat = read.table(gt.file, sep='\t', stringsAsFactors=F, header=T)
  colnames(gt.mat) = gsub('.*\\.run1\\.(.*)_12h.*', '\\1', colnames(gt.mat))
  id.num = paste(gt.mat[, 'CHROM'], gt.mat[, 'POS'], sep='.')
  rownames(gt.mat) = id.num
  gt.mat = gt.mat[, setdiff(colnames(gt.mat), c('CHROM', 'POS'))]

  #sub 00, 01, 11 -> 0,1,2
  gt.num.mat = as.data.frame(lapply(gt.mat, function(ind.gt){ind.gt = gsub('0/0', '0', ind.gt); ind.gt = gsub('0/1', '1', ind.gt); ind.gt = gsub('1/1', '2', ind.gt); ind.gt = as.integer(ind.gt); return(ind.gt);}))
  rownames(gt.num.mat) = rownames(gt.mat)
  colnames(gt.num.mat) = gsub('^X', '', colnames(gt.num.mat))
  gt.mat = gt.num.mat

  return(gt.mat)
}

read.dp <- function(dp.file){
  dp.mat = read.table(dp.file, sep='\t', stringsAsFactors=F, header=T)
  colnames(dp.mat) = gsub('.*\\.run1\\.(.*)_12h.*', '\\1', colnames(dp.mat))
  colnames(dp.mat) = gsub('^X', '', colnames(dp.mat))
  id.num = paste(dp.mat[, 'CHROM'], dp.mat[, 'POS'], sep='.')
  rownames(dp.mat) = id.num
  dp.mat = dp.mat[, setdiff(colnames(dp.mat), c('CHROM', 'POS'))]

  return(dp.mat)
}

#Execute
main(dp.file, gt.file, dp.min)
