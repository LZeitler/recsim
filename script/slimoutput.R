library(dplyr)
library(ggplot2)
library(cowplot)

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

df <- getoutput('out-equi.txt')
head(df)

ggplot(df)+
    geom_line(aes(gen,neutral),color='black')+
    geom_line(aes(gen,del),color='red')+
    geom_line(aes(gen,ben),color='green')+
    labs(x='Generation',y='Segregating sites')

## ideas: density in regions for mutations classes (count/MBp) after burn in
## show allele freq for freq type and positions qith rec background (implement in equi.E)
