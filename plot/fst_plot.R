source('~/code/r/source_me.R')
source('~/pro/recsim/plot/funs.R')
setwd('~/pro/recsim/300_analyses/')

files <- list.files(path='/mnt/evo_euler/data/300/fst/',pattern="*.txt")
genes <- fread('annotations/yant_genes.csv',data.table=F) %>% filter(meiosis==T)


file <- files[1]
dt <- fread(file,data.table=F)

