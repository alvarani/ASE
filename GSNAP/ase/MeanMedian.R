


##mean & median for GSNAP and degner
#Gsnap:
load("/proj/b2012046/rani/analysis/gsnap/ase/ase.res.RData") 
#degner:
load("/proj/b2012046/rani/analysis/degner/ase/ase.res.RData")
ls(ase.res)
sapply(ase.res[], mean, c(1))# for all mean for all columns selected alt.frac and then calculated the average manually
#median
sapply(ase.res[[1]], median, c(1)) # so as for 2 till 16 and taken alt.frac and then taken the average for meadian
