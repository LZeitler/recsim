library(dplyr)
library(data.table)

run <- '142017241019'
rundir <- paste0('/cluster/scratch/zeitlerl/slim/',run)

setwd(rundir)

par <- fread('parspace.txt',data.table=F)


fnames <- system("find . -name 'qtl3-phenotypes.txt'",intern=T)

dir.create('all_pheno',recursive=T)

for (r in par[,'parcomb']) dir.create(paste0('all_pheno/',r)) # create new dirs

for (r in fnames){                                 # move files

    parcomb <- substr(r,3,6)                       # parcomb number
    seed <- substr(r,8,20)                         # seed
    basen <- 'qtl3-phenotypes'
    
    system(
        paste0('cp ',r,' all_pheno/',parcomb,'/',basen,'_',seed,'.txt')
    )
    
}

setwd('all_pheno')

big <- data.frame()
for (r in list.files()){
    s <- lapply(list.files(r), function(x) data.frame(fread(paste0(r,'/',x),data.table=F),par=r,seed=substr(x,17,29)))
    big <- rbind(big,do.call(rbind,s))
}

fwrite(big,paste0('~/pro/recsim/output/pheno-',run,'.txt'))
