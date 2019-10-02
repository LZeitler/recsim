makemap <-                              # function to generate a map with 3 seperate rec regions
    function(segments=c(.4,.6),length=100000,means=c(5e-8,1e-9,5e-8),sds=c(1e-8,1e-8,1e-8),sites=1000){
        rellen <- c(sites*segments[1],
                    sites*segments[2]-sites*segments[1],
                    sites-sites*segments[2])
        t <- data.frame(rate=abs(c(rnorm(rellen[1],means[1],sds[1]),
                                   rnorm(rellen[2],means[2],sds[2]),
                                   rnorm(rellen[3],means[3],sds[3]))),
                        pos= round(c(seq(0,floor(length*segments[1]),length.out=rellen[1]),
                                     seq(ceiling(length*segments[1]),floor(length*segments[2]),length.out=rellen[2]),
                                     seq(ceiling(length*segments[2]),length-1,length.out=rellen[3]))))
        ## rellen
        t
    }

source('~/pro/300_analyses/script/r/recomb.funs.R')
source('~/code/r/source_me.R')

set.seed(1345)
rmap <- makemap(segments = c(.2,.8), means = c(5e-8,5e-10,5e-8), sds = c(2e-8,2e-8,2e-8))
a <- qplot(rmap$pos,rmap$rate)+geom_smooth()
rmap$crate <- mapply(function(r,w) haldane(rec=r,wi=w), rmap$rate,rmap$pos)
rmap$ccrate <- cumsum(rmap$crate)
b <- qplot(rmap$pos,rmap$ccrate)
plot_grid(a,b,nrow=2)
max(rmap$ccrate)                        # 150

set.seed(11222345)
rmap <- makemap(segments = c(.4,.6), means = c(3.2e-8,3e-10,3.2e-8), sds = c(2e-8,2e-8,2e-8))
a <- qplot(rmap$pos,rmap$rate)+geom_smooth()
rmap$crate <- mapply(function(r,w) haldane(rec=r,wi=w), rmap$rate,rmap$pos)
rmap$ccrate <- cumsum(rmap$crate)
b <- qplot(rmap$pos,rmap$ccrate)
plot_grid(a,b,nrow=2)
max(rmap$ccrate)                        # 150


############### with less variation
set.seed(12345)
rmap <- makemap(segments = c(.2,.8), means = c(6.4e-8,5e-10,6.4e-8), sds = c(1e-8,1e-8,1e-8))
a <- qplot(rmap$pos,rmap$rate)+geom_smooth()
rmap$crate <- mapply(function(r,w) haldane(rec=r,wi=w), rmap$rate,rmap$pos)
rmap$ccrate <- cumsum(rmap$crate)
b <- qplot(rmap$pos,rmap$ccrate)
plot_grid(a,b,nrow=2)
max(rmap$ccrate)                        # 150

fwrite(select(rmap,pos,rate),'/home/leo/projects/recsim/slim/maps/150cm-2-8-lowvar.bed',quote=F,sep='\t',col.names=F)

set.seed(1235)
rmap <- makemap(segments = c(.4,.6), means = c(3.5e-8,1e-10,3.5e-8), sds = c(1e-8,1e-8,1e-8))
a <- qplot(rmap$pos,rmap$rate)+geom_smooth()
rmap$crate <- mapply(function(r,w) haldane(rec=r,wi=w), rmap$rate,rmap$pos)
rmap$ccrate <- cumsum(rmap$crate)
b <- qplot(rmap$pos,rmap$ccrate)
plot_grid(a,b,nrow=2)
max(rmap$ccrate)                        # 150

fwrite(select(rmap,pos,rate),'/home/leo/projects/recsim/slim/maps/150cm-4-6-lowvar.bed',quote=F,sep='\t',col.names=F)


############### with less variation, a bigger chromosome

easymap <- function(r,rcenr,segments=c(.2,.6,.2),wi=1000,chromlen=1e7){
    t <- data.frame(rec=c(rep(r,segments[1]*wi),rep(rcenr,segments[2]*wi),rep(r,segments[3]*wi)),
                    wi=seq(chromlen/wi,chromlen,length.out=wi))
    t$cm <- c(rep(haldane(r,chromlen/wi),segments[1]*wi),
              rep(haldane(rcenr,chromlen/wi),segments[2]*wi),
              rep(haldane(r,chromlen/wi),segments[3]*wi))
    t$ccm <- cumsum(t$cm)
    cat(paste('Chromosome length is',round(max(t$ccm),2),'cM \n'))
    return(t)
}

easymap2 <- function(r,rcenr,segments=c(.2,.6,.2),wi=1000,chromlen=1e7,sd=0,sdcenr=0){
    t <- data.frame(rec=c(rnorm(segments[1]*wi,r,sd),rnorm(segments[2]*wi,rcenr,sdcenr),rnorm(segments[3]*wi,r,sd)),
                    wi=seq(chromlen/wi,chromlen,length.out=wi))
    t$rec[t$rec<0] <- 0                # set negative recombination to 0
    t$cm <- mapply(function(re,w) haldane(rec=re,wi=w), t$rec, chromlen/wi)
    t$ccm <- cumsum(t$cm)
    cat(paste('Chromosome length is',round(max(t$ccm),2),'cM \n'))
    return(t)
}

t <- easymap(2e-7,1e-8,c(.4,.2,.4))
qplot(t$wi,t$ccm)
qplot(t$wi,t$rec)
qplot(t$wi,t$cm)

## chromsome lengths approx 150 cM

set.seed(123)
a <- easymap2(1.9e-7,1e-8,c(.4,.2,.4),sd=1e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/150cm-4-2-easy-lovar.bed',quote=F,sep='\t',col.names=F)


set.seed(1234)
a <- easymap2(2.5e-7,1e-8,c(.3,.4,.3),sd=1e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/150cm-3-4-easy-lovar.bed',quote=F,sep='\t',col.names=F)


set.seed(1234)
a <- easymap2(3.6e-7,1e-8,c(.2,.6,.2),sd=1e-8,sdcenr=3e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/150cm-2-6-easy-lovar.bed',quote=F,sep='\t',col.names=F)



set.seed(12345)
rmap <- makemap(segments = c(.2,.8), means = c(6.4e-8,5e-10,6.4e-8), sds = c(1e-8,1e-8,1e-8),length=1E6)
a <- qplot(rmap$pos,rmap$rate)+geom_smooth()
rmap$crate <- mapply(function(r,w) haldane(rec=r,wi=w), rmap$rate,rmap$pos)
rmap$ccrate <- cumsum(rmap$crate)
b <- qplot(rmap$pos,rmap$ccrate)
plot_grid(a,b,nrow=2)
max(rmap$ccrate)                        # 150

fwrite(select(rmap,pos,rate),'/home/leo/projects/recsim/slim/maps/150cm-2-8-lowvar.bed',quote=F,sep='\t',col.names=F)

set.seed(1235)
rmap <- makemap(segments = c(.4,.6), means = c(3.5e-8,1e-10,3.5e-8), sds = c(1e-8,1e-8,1e-8))
a <- qplot(rmap$pos,rmap$rate)+geom_smooth()
rmap$crate <- mapply(function(r,w) haldane(rec=r,wi=w), rmap$rate,rmap$pos)
rmap$ccrate <- cumsum(rmap$crate)
b <- qplot(rmap$pos,rmap$ccrate)
plot_grid(a,b,nrow=2)
max(rmap$ccrate)                        # 150

fwrite(select(rmap,pos,rate),'/home/leo/projects/recsim/slim/maps/150cm-4-6-lowvar.bed',quote=F,sep='\t',col.names=F)

## chromsome lengths approx 75 cM

set.seed(123)
a <- easymap2(9e-8,1e-8,c(.4,.2,.4),sd=1e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/75cm-4-2-easy-lovar.bed',quote=F,sep='\t',col.names=F)


set.seed(123)
a <- easymap2(1.18e-7,1e-8,c(.3,.4,.3),sd=1e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/75cm-3-4-easy-lovar.bed',quote=F,sep='\t',col.names=F)


set.seed(1234)
a <- easymap2(1.62e-7,1e-8,c(.2,.6,.2),sd=1e-8,sdcenr=3e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/75cm-2-6-easy-lovar.bed',quote=F,sep='\t',col.names=F)


## chromsome lengths approx 100 cM with more realistic values

set.seed(123)
a <- easymap2(1.12e-7,5e-8,c(.4,.2,.4),sd=2e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-4-2-easy-realv.bed',quote=F,sep='\t',col.names=F)


set.seed(1224)
a <- easymap2(1.33e-7,5e-8,c(.3,.4,.3),sd=2e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-3-4-easy-realv.bed',quote=F,sep='\t',col.names=F)


set.seed(1236)
a <- easymap2(1.75e-7,5e-8,c(.2,.6,.2),sd=2e-8,sdcenr=1e-8,wi=1000)
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-2-6-easy-realv.bed',quote=F,sep='\t',col.names=F)


## chromsome lengths approx 100 cM without any variation

set.seed(123)
a <- easymap2(1.12e-7,5e-8,c(.4,.2,.4))
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-4-2-easy-novar.bed',quote=F,sep='\t',col.names=F)


set.seed(1224)
a <- easymap2(1.33e-7,5e-8,c(.3,.4,.3))
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-3-4-easy-novar.bed',quote=F,sep='\t',col.names=F)


set.seed(1236)
a <- easymap2(1.75e-7,5e-8,c(.2,.6,.2))
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-2-6-easy-novar.bed',quote=F,sep='\t',col.names=F)


## chromsome lengths approx 100 cM without any variation and no recombination in center

set.seed(123)
a <- easymap2(1.25e-7,0,c(.4,.2,.4))
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-4-2-easy-novar-nocent.bed',quote=F,sep='\t',col.names=F)


set.seed(1224)
a <- easymap2(1.67e-7,0,c(.3,.4,.3))
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-3-4-easy-novar-nocent.bed',quote=F,sep='\t',col.names=F)


set.seed(1236)
a <- easymap2(2.5e-7,0,c(.2,.6,.2))
qplot(a$wi,a$ccm)
qplot(a$wi,a$rec)

fwrite(select(a,wi,rec),'/home/leo/projects/recsim/slim/maps/100cm-2-6-easy-novar-nocent.bed',quote=F,sep='\t',col.names=F)

####################
## reverse
sfun <- function(rec) 100*-0.5*log(1-2*rec*1000)
ss <- seq(0,0.0001,0.000000001)
ssd <- data.frame(ss) %>% mutate(sfun=sfun(ss))
ssd[which.min(abs((tlen/len*sites)-ssd$sfun)),1] # per bp recombination rate

