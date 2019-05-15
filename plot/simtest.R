source('~/code/r/source_me.R')
source('funs.R')

dt <- fread('../script/simtest_out.txt',data.table=F)
a <- ggplot(dt,aes(gen,freq))+
    geom_point()+
    geom_line()+
    coord_cartesian(ylim=c(0,1))
a

dt <- fread('../script/simtest_out_rep.txt',data.table=F)
## dt <- bind_rows(guessstart(dt),dt)
dt <- bind_rows(knowstart('../script/simtest_stats.txt'),dt)
head(filter(dt,rep==1))
a <- ggplot(dt,aes(gen,freq,group=factor(rep),color=factor(rep)))+
    ## geom_point()+
    geom_line()+
    coord_cartesian(ylim=c(0,1))
a
