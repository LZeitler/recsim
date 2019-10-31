## QTL4 submitter
## submits array jobs with desired parameters as job array
## run from login shell

library(dplyr)
library(data.table)

mydate <- function(){
    paste0(format(Sys.time(),format="%H%M%S"),format(Sys.Date(),format="%d%m%y"))
}

wdir <- paste0('/cluster/scratch/zeitlerl/slim/',mydate())
dir.create(wdir,recursive=T)
setwd(wdir)

slimpath <- "~/programs/SLiM_build/slim"
scriptpath <- "~/pro/recsim/script/eidos/qtl4-all.E"

mparams <- function(params){
    m <- lapply(params,function(x) length(unique(x)))
    names(params)[which(m>1)]
}
paramsspace <- function(params){
    mpar <- mparams(params)
    if (length(mpar)>1){
        l <- sapply(mpar,function(x) eval(parse(text=x)))
        data.frame(expand.grid(l, stringsAsFactors = F),select(params,-mpar), stringsAsFactors = F)
    } else {
        data.frame(params,stringsAsFactors = F)
    }
    
}
slimrunner <- function(seed, script, define, dryrun=F){
    if (dryrun){
        paste(slimpath,
              ## '-s',seed,
              define,script)
    } else {
        system2(slimpath, args = c(
                              ## "-s", seed,
                              define,
                              shQuote(script)), stdout=T, stderr=T)
    }
}
runner <- function(workdir,run){
    path <- paste0(workdir,'/',run)
    dir.create(path)
    setwd(path)
}

## -d[efine] <def>  : define an Eidos constant, such as "mu=1e-7"


## define replicates
reps <- 100


## define parameters
mdel=c(-0.005,-0.01,-0.05)
mben=c(0.05)
nneutral=c(1)
ndel=c(0.05)
nben=c(0.01)                                                       
n=c(50,100,500)
opt=c(10)                                                          
sigma=c(5)                                     
scale=dnorm(0.0, 0.0, sigma)
mu=1e-7
hdel=c(0,.5)
hben=c(.5)
mapfile=c('\u5c\u5c\u27~/pro/recsim/slim/maps/050cm-2-6-easy-016.bed\u5c\u5c\u27',
          '\u5c\u5c\u27~/pro/recsim/slim/maps/050cm-3-4-easy-016.bed\u5c\u5c\u27',
          '\u5c\u5c\u27~/pro/recsim/slim/maps/050cm-4-2-easy-016.bed\u5c\u5c\u27')
          

################
## run everything
params <- list(
    mdel=mdel,
    mben=mben,
    nneutral=nneutral,
    ndel=ndel,
    nben=nben,
    n=n,
    opt=opt,
    sigma=sigma,
    scale=scale,
    mu=mu,
    hdel=hdel,
    hben=hben,
    mapfile=mapfile
)

## parspace <- paramsspace(params)
parspace <- expand.grid(params)
cat('Calculated',nrow(parspace),'parameter combinations\n')
parspace <- parspace %>% mutate(parcomb=1000+1:nrow(parspace))

fwrite(parspace,paste0('parspace.txt'))

cat('Written parameter space table.\n')
cat(wdir,'\n\n')

## out <- as.vector(NULL)
for (r in 1:nrow(parspace)){

    comb <- parspace[r,]
    combn <- parspace$parcomb[r]
    ## for (i in 1E5:(1E5+reps)){
    defstr <- paste(sapply(1:ncol(comb),function(x)
        paste0('-d ',names(comb)[x],'=',comb[,x])),collapse=' ')
    ## run <- paste0(combn,i)
    run <- combn
    slimcmd <- slimrunner(run,scriptpath,defstr,T)
    bsubcmd <- paste(c("bsub -J 'qtl4",run,"[1-",reps,
                       "]%50' -n 1 -W 12:00 -R 'rusage[mem=4000]' -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr "),
                     collapse = '')

    ## replicates are submitted as array jobs, bsub cmds are written to file and then submitted
    bsubs <- paste0(bsubcmd,shQuote(slimcmd,'cmd')) # map without \

    fwrite(list(bsubs),paste0(wdir,'/','bsubcmds.sh'),quote=F,col.names=F,sep='\t',append=T)
    cat('slim-runner.R wrote commands to file bsubcmds.sh.\n')
    
    ## functional command looks like this
    ## bsub -J 'qtl41001[1-1]' -n 1 -W 12:00 -R 'rusage[mem=4000]' -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr "~/programs/SLiM_build/slim -d mdel=-0.001 -d mben=0.01 -d nneutral=1 -d ndel=0.1 -d nben=0.01 -d n=10 -d opt=10 -d sigma=5 -d scale=0.0797884560802865 -d mu=1e-07 -d mapfile=\\'~/pro/recsim/slim/maps/100cm-4-2-easy-novar-nocent.bed\\' -d parcomb=1001 ~/pro/recsim/script/eidos/qtl4-all.E"
    
}

system('bash bsubcmds.sh')
cat('slim-runner.R submitted all jobs.\n')

## accumulator <- integer(reps)
## for (iter in 1:reps){
## 	# collect output from running the script once
## 	output <- slimrunner(iter, scriptpath)
	
## 	# find a pattern line we know is in the output just before the data we want
## 	leadLineIndex <- which(grepl("^Fraction with p1 ancestry, by position:$", output))
	
## 	# extract the actual numeric data from the data line's string
## 	valuesLine <- output[leadLineIndex + 1]
## 	valuesList <- strsplit(valuesLine, " ", fixed=T)
## 	valuesString <- unlist(valuesList)
## 	values <- as.numeric(valuesString)
	
## 	# add the data into our accumulator buffer
## 	accumulator <- accumulator + values
## }
