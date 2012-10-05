#Reads annotated vars, adds column names and dumps as rdata object
#
#USAGE: Rscript annovar_flatfile2rdata.R varcalls.infile varcalls.outfile

argv = commandArgs(trailingOnly=TRUE)
varcalls.infile = argv[1]
varcalls.outfile = argv[2]

#varcalls.infile = '/proj/b2011075/analysis/hs_pipe/outdata/run1/ase/annovar/annovar_all_variants.txt'
#varcalls.outfile = '/proj/b2011075/analysis/hs_pipe/outdata/run1/ase/annovar/annovar_all_variants.wheader.txt'
#gene rel and gene annot read from annovar files: *.variant_function. I think annovar annot by ref_gene db.
first.colnames = c('chr', 'start', 'stop', 'refvar', 'observed var', 'gene relation', 'gene annot') #7
fixed.dbs = c('var type', 'aa change', 'ConservedReg', 'SegDup', '1000g2010nov', 'Dbsnp_1', 'Dbsnp_2') #7
last.dbs = c('cg46', 'Avsift', 'PP2', 'dgv', 'mirna', 'tfbs', 'ljb_mt', 'ljb_lrt', 'gwascatalog', 'mirnatarget', 'omimgene', '1000g2011may', '1000g2010nov.EUR.snps', '1000g2010nov.EUR.indels') #14
#FLAT file:
#6: gene rel
#7: gene annot
#8: empty
#9: var type
main <- function(varcalls.infile, varcalls.outfile, first.colnames, fixed.dbs, last.dbs){

  #Read varcalls
  varcalls = read.varcalls(varcalls.infile, first.colnames, fixed.dbs, last.dbs)
  n.vars = nrow(varcalls)
  print(sprintf('n vars: %i', n.vars))

  #Coerce to numeric cols so that get.lt.ind and get.gt.ind works!
  char2num.cols = c(colnames(varcalls)[grep('^1000g', colnames(varcalls))], 'cg46', 'Avsift', 'PP2', 'ljb_mt', 'ljb_lrt', colnames(varcalls)[grep('\\.dp$', colnames(varcalls))])
  for(jcol in char2num.cols){
    varcalls[, jcol] = as.numeric(varcalls[, jcol])
  }

  #add indel column
  varcalls = add.indel.col(varcalls)  
  print(colnames(varcalls))

  #add rownames: chr.start.obsvar
  site = paste(varcalls[, 'chr'], varcalls[, 'start'], varcalls[, 'observed var'], sep = '.')
  rownames(varcalls) = site
  
  #Write new table with added depths
  write.cols = setdiff(colnames(varcalls), 'chr.start')
  varcalls = varcalls[, write.cols]
  write.table(varcalls, file = varcalls.outfile, sep='\t', row.names=F, quote=F)

  #dump as RData object
  varcalls.rdata.file = paste(varcalls.outfile, 'RData', sep = '.')
  save(varcalls, file = varcalls.rdata.file)
}

read.varcalls <- function(varcall.file, first.colnames, fixed.dbs, last.dbs){

  varcalls = read.table(varcall.file, sep='\t', stringsAsFactors = FALSE, quote='')
  #rm last col (empty)
  varcalls = varcalls[, 1:(ncol(varcalls) - 1)]
  
  last.colnames = c(fixed.dbs, last.dbs)  
  colnames(varcalls)[1:length(first.colnames)] = first.colnames
  colnames(varcalls)[(ncol(varcalls) - length(last.colnames) + 1):ncol(varcalls)] = last.colnames

  varcalls = varcalls[, c(first.colnames, last.colnames)] #order cols

  #Rm 'chr' prefix from chromosome name
  varcalls[, 'chr'] = gsub('chr', '', varcalls[, 'chr'])

  #Rm db prefix from 1000g
  varcalls[, '1000g2010nov'] = gsub('.*;', '', varcalls[, '1000g2010nov'])
  
  #Rm db prefix
  for(db in last.dbs){
    varcalls[, db] = gsub('.*;', '', varcalls[, db])
  }

  #Rm 'Name' prefix
  for(db in last.dbs){
    varcalls[, db] = gsub('Name=', '', varcalls[, db])
  }
    
  return(varcalls)
}

add.indel.col <- function(vars){
#NB: assumes that if not SNV than it has to be an indel

  #get indel indices
  n.vars = nrow(vars)
  snv.ind = get.snv.ind(vars)
  indel.ind = setdiff(1:n.vars, snv.ind)

  #add indel col
  indel = rep(0, n.vars)
  indel[indel.ind] = 1  
  vars = cbind(vars, indel)
  
  return(vars)
}

get.snv.ind <- function(vars){

  obsvars = vars[, 'observed var']
  refvars = vars[, 'refvar']
  bases = c('A', 'C', 'G', 'T')

  #get non-insertions
  nonins.ind = integer(0)
  for(jbase in bases){
    nonins.ind = unique(c(nonins.ind, which(obsvars == jbase)))
  }

  #get non-deletions
  nondel.ind = integer(0)
  for(jbase in bases){
    nondel.ind = unique(c(nondel.ind, which(refvars == jbase)))
  }

  #get snvs
  snv.ind = intersect(nonins.ind, nondel.ind)

  return(snv.ind)
}

#Execute
main(varcalls.infile, varcalls.outfile, first.colnames, fixed.dbs, last.dbs)
