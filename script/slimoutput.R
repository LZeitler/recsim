library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)

getoutput <- function(file){
    system(paste0('grep Fixed: ../slim/', file, ' > ../slim/', file,'-n.txt'))
    xx <- data.table::fread(paste0('../slim/', file,'-n.txt'),data.table=F) %>% t() %>% as.vector()
    xx <- as.numeric(xx[2:length(xx)])
    xx
}

getoutput <- function(file){
    system(paste0('grep "Mean phenotype:" ../slim/', file, ' > ../slim/', file,'-n.txt'))
    xx <- data.table::fread(paste0('../slim/', file,'-n.txt'),data.table=F) %>% t() %>% as.vector()
    xx <- as.numeric(xx[2:length(xx)])
    xx
}



yy <- getoutput('out6.txt')

yy[yy!=0]
length(yy[yy<0])

qplot(yy,bins=10)
qplot(round(yy[yy!=0],3),bins=12)
(yp <- qplot(round(yy[yy<0],3),bins=12))


xx <- getoutput('out5.txt')

qplot(round(xx[xx!=0],3),bins=12)
(xp <- qplot(round(xx[xx<0],3),bins=12))

xx <- getoutput('qtl1.txt')
qplot(1:length(xx),xx)

plot_grid(yp,xp,ncol=1)




####################

## new script with fixed mutation types (burn-in)
getoutput <- function(file){
    system(paste0("awk '/gen\t/ {seen = 1} seen {print}' ", file, " > ", file, "-output.txt"))
    xx <- data.table::fread(paste0(file,'-output.txt'),data.table=F)
    xx
}

setwd('../slim')

df <- getoutput('out-qtl3-equi.txt')
head(df)

ggplot(df)+
    geom_line(aes(gen,neutral),color='black')+
    geom_line(aes(gen,del),color='red')+
    geom_line(aes(gen,ben),color='green')+
    labs(x='Generation',y='Segregating sites')

## site pi in vcf
pi <- fread('samplevcf-equi-sitepi.sites.pi',data.table=F)
sel <-fread('samplevcf-equi-selcoeff.INFO',data.table=F)
go <-fread('samplevcf-equi.INFO',data.table=F)
pisel <- inner_join(pi,sel) %>%
    mutate(mtype=ifelse(round(S,3)!=0,
                 ifelse(round(S,3)>0,'pos','neg'),'0'),
           region=ifelse(POS<3E6|POS>7E6,'high','low'))
pisel <- inner_join(pisel,go)

ggplot(pisel)+
    geom_point(aes(POS,PI,color=S>0))

hist(filter(pisel,S!=0)$S)

pisel$window <- cut(pisel$POS,seq(0,1E7,10000))
pisel <- pisel %>% group_by(window) %>%
    mutate(meanpi=mean(PI),
           nmut=length(PI)) %>%
    group_by(window,mtype) %>%
    mutate(nmuttype=length(PI)) %>% data.frame

ggplot(pisel)+
    geom_point(aes(POS,meanpi,color=S>0))

ggplot(pisel)+
    geom_jitter(aes(POS,nmut,color=mtype))+
    geom_vline(aes(xintercept=3E6))+
    geom_vline(aes(xintercept=7E6))

ggplot(pisel)+
    geom_boxplot(aes(region,PI,color=mtype))

ggplot(pisel)+
    geom_boxplot(aes(region,nmuttype,color=mtype))

ggplot(pisel)+
    geom_boxplot(aes(region,10000-GO,color=mtype))

## ideas: density of mutations in regions for mutations classes (count/MBp) after burn in
## show allele freq for freq type and positions qith rec background (implement in equi.E)


######################
######################
## QTL with maps and adaptive walk

walk <- fread('1871742102694/qtl3-phenotypes.txt',data.table=F)

ggplot(walk)+
    geom_line(aes(generation,mphenotype))

walk3 <- fread('1847810905742/qtl3-phenotypes.txt',data.table=F)  # burn in with d b n with selection effect
walk2 <- fread('1854323474551/qtl3-phenotypes.txt',data.table=F) # burn in with d b n w/o selection effect
walk4 <- fread('1856494609495/qtl3-phenotypes.txt',data.table=F) # neutral burn in
walk1 <- fread('1870891446584/qtl3-phenotypes.txt',data.table=F) # no burn in

walkes <- rbind(data.frame(walk3, sim='with n/b/d standing var'),
                data.frame(walk4, sim='with n standing var'),
                data.frame(walk2, sim='with n/b/d standing var, no S'),
                data.frame(walk1, sim='w/o standing var'))

ggplot(walkes)+
    geom_line(aes(generation,mphenotype,color=sim))

eff <- fread('1750863695504/qtl2-effects.txt',data.table=F)

ggplot(eff)+
    geom_histogram(aes(effect),binwidth = .05)

ggplot(eff)+
    geom_histogram(aes(freq,fill=effect>0),position='dodge',binwidth=.05)
