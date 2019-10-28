library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=T)
if (length(args)==0){
    stop('No run number supplied. Please add run number as argument to Rscript command. No default\n
    	  Usage: Rscript slim-collector.R --args <run number> <which output>\n\n')
} else {
    run <- as.character(args[2])
    cat('Run: ',run,'\n\n')
    this <- as.character(args[3])
    cat('Selected output',this,'\n\n')
}

rundir <- paste0('/cluster/scratch/zeitlerl/slim/',run)

setwd(rundir)

par <- fread('parspace.txt',data.table=F)

if (this=='qtl3-effects.txt') odir <- 'all_beneffects'
if (this=='qtl3-deleffects.txt') odir <- 'all_deleffects'
if (this=='qtl3-phenotypes.txt') odir <- 'all_pheno'

fnames <- system(paste0("find . -name ", this),intern=T)

dir.create(odir,recursive=T)

for (r in par[,'parcomb']) dir.create(paste0(odir,'/',r)) # create new dirs

for (r in fnames){                                 # move files

    parcomb <- substr(r,3,6)                       # parcomb number
    seed <- substr(r,8,20)                         # seed
    basen <- gsub('.txt','',this)
    
    system(
        paste0('cp ',r,' ',odir,'/',parcomb,'/',basen,'_',seed,'.txt')
    )
    
}

setwd(odir)

big <- data.frame()
for (r in list.files()){
    s <- lapply(list.files(r), function(x) data.frame(fread(paste0(r,'/',x),data.table=F),par=r,seed=regmatches(r,regexec("_(.*?).txt",r))[[1]][2]))
    big <- rbind(big,do.call(rbind,s))
}

fwrite(big,paste0('~/pro/recsim/output/',gsub('all_','',odir),'-',run,'.txt'))
cat('Wrote output.\n')

system(paste0('cp ', rundir, '/parspace.txt ~/pro/recsim/output/parspace-', run, '.txt'))
cat('Copied parameter space.\n')
