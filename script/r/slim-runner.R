## QTL4 submitter
## submits array jobs with desired parameters as job array
## run from login shell

mydate <- function(){
    paste0(format(Sys.time(),format="%H%M%S"),format(Sys.Date(),format="%d%m%y"))
}

wdir <- paste0('$SCRATCH/slim/',mydate())
system2(paste0('mkdir -p ',wdir)
setwd(wdir)

slimpath <- "~/programs/SLiM_build/slim"
scriptpath <- "~/pro/recsim/script/eidos/qtl4-all.E"

mparams <- function(params){
    m <- sapply(1:ncol(params),function(x) nrow(unique(params[x])))
    names(params)[which(m>1)]
}
paramsspace <- function(params){
    mpar <- mparams(params)
    l <- sapply(mpar,function(x) eval(parse(text=x)))
    data.frame(expand.grid(l),select(params,-mpar))
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
    mkdir(path)
    setwd(path)
}

## -d[efine] <def>  : define an Eidos constant, such as "mu=1e-7"


## define replicates
reps <- 1


## define parameters
mdel=c(-0.001)
mben=c(0.01)
nneutral=c(1)
ndel=c(0.1)
nben=c(0.01)                                                       
n=c(10,100,1000)                                                        
opt=c(10)                                                          
sigma=c(5)                                     
scale=dnorm(0.0, 0.0, sigma)
mu=1e-7
mapfile=c("~/pro/recsim/slim/maps/100cm-4-2-easy-novar-nocent.bed",
          "~/pro/recsim/slim/maps/100cm-3-4-easy-novar-nocent.bed",
          "~/pro/recsim/slim/maps/100cm-2-6-easy-novar-nocent.bed",
          "~/pro/recsim/slim/maps/100cm-4-2-easy-novar.bed",
          "~/pro/recsim/slim/maps/100cm-3-4-easy-novar.bed",
          "~/pro/recsim/slim/maps/100cm-2-6-easy-novar.bed"
          )


################
## run everthing
params <- data.frame(
    mdel,
    mben,
    nneutral,
    ndel,
    nben,
    n,
    opt,
    sigma,
    scale,
    mu,
    mapfile
)

parspace <- paramsspace(params)
parspace.out <- parspace %>% mutate(parcomb=1000+1:nrow(parspace))

fwrite(parspace.out,paste0('parspace-run.txt')))

## out <- as.vector(NULL)
for (r in 1:nrow(parspace)){
    comb <- parspace[r,]
    combn <- parspace.out$parcomb[r]
    ## for (i in 1E5:(1E5+reps)){
    defstr <- paste(sapply(1:ncol(comb),function(x)
        paste0('-d ',names(comb)[x],'=',comb[,x])),collapse=' ')
    ## run <- paste0(combn,i)
    run <- combn
    slimcmd <- slimrunner(run,scriptpath,defstr,T)
    runner(wdir,run)
    ## replicates are submitted as array jobs
    system2(noquote(paste0('bsub -J "qtl4',run,'[1-',reps,
                           ']" -n 1 -W 12:00 -R "rusage[mem=4000]" -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr "',
                           slimcmd,'"')))
    ## out <- c(out,slimcmd)
    ## }
}

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
