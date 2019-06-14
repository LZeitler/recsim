source('~/code/r/source_me.R')
source('~/pro/recsim/plot/funs.R')
setwd('~/pro/300_analyses/')

nfs <- '/mnt/evo_euler/data/300/fst/'                 # local
## nfs <- '/nfs/nas22/fs2201/biol_impb_group_evo_euler/' # remote

files <- list.files(path=nfs,pattern="*.txt")
genes <- fread('annotations/yant_genes.csv',data.table=F) %>% filter(meiosis=='T') %>%
    rename(chr=Chr,ATG_ID=`ATG ID`) %>% select(c(1:3,5)) %>% mutate(ATG_ID=make.unique(ATG_ID))
## genes.m <- melt(genes,c('chr','ATG ID'),value.name = 'pos')

fst.plot <- function(data,...){
    p <- ggplot(data,aes(pos,FstSNP))+
        ## geom_point()+
        stat_summary_bin(binwidth = 100000,
                         fun.y = mean,
                         geom = 'line',
                         size = .5)+
        facet_wrap(.~chr,nrow=1,scales='free_x')+
        coord_cartesian(ylim = c(0,1))
    p
}

file <- files[2]
dt <- fread(paste0(nfs,file),data.table=F)


dt <- filter(dt,chr%in%c(1,2,4,6))

## dt <- top_n(dt,round(nrow(dt)/10))

fst.plot(dt)

fst.plot(dt)+
    geom_vline(data = genes,aes(xintercept=start),color='green',alpha=.5)+
    geom_vline(data = genes,aes(xintercept=end),color='red',alpha=.5)

## make windows for FST
start <- dt%>%group_by(chr)%>%summarise(start=min(pos),end=max(pos))
win <- sapply(1:nrow(start), function(x) seq(start$start[x],start$end[x],10000))
wi <- data.frame()
for (i in 1:nrow(start)) {wi <- rbind(data.frame(chr=i,pos=win[[i]],window=win[[i]]),wi)}
dtw <- full_join(dt,wi) %>% arrange(chr,pos) %>% fill(window) %>% data.frame 
fstw <- dtw %>% group_by(chr,window) %>% summarise(FstSNP=mean(FstSNP,na.rm=T),nsnp=n()) %>%
    arrange(chr,window) %>% data.frame
summary(fstw$nsnp)
summary(fstw$FstSNP)

## filter out bad windows
fstw <- fstw %>% filter(FstSNP>0, !is.na(FstSNP), nsnp>=5)

## plot it 
fst.plot(dtw)+
    geom_vline(data = genes,aes(xintercept=start),color='green',alpha=.5)+
    geom_vline(data = genes,aes(xintercept=end),color='red',alpha=.5)


## find FST for genes with a random locations (alternative: circular shuffle)
setDT(dt);setDT(genes)

dtg <- genes[dt,on= .(chr==chr,         # SNPs within genes
                      start<=pos,
                      end>=pos),
          nomatch=0,
          .(chr,pos,ATG_ID,FstSNP)]

randpos <- lapply(1:nrow(genes),function(x) data.frame(ATG_ID=genes[x,4], end=genes[x,3]-genes[x,2],
                                                       sample_n(filter(dt,chr==unlist(genes[x,1],use.names = F),
                                                                       !is.na(FstSNP)),1000)))
randpos <- do.call(rbind,randpos) %>% mutate(start=pos,end=start+end) %>% select(chr,start,end,ATG_ID) %>% setDT

dtr <- randpos[dt,on= .(chr==chr,         # SNPs with random pos but length of original genes
                      start<=pos,
                      end>=pos),
          nomatch=0,
          .(chr,pos,ATG_ID,FstSNP)]

(av.genes <- dtg%>%group_by(ATG_ID)%>%summarise(meanFST.genic=mean(FstSNP,na.rm=T)))
(av.random <- dtr%>%group_by(ATG_ID)%>%summarise(meanFST.random=mean(FstSNP,na.rm=T)))
(av.genome <- dt%>%group_by(chr)%>%summarise(meanFST.genome=mean(FstSNP,na.rm=T)))
out <- inner_join(genes,av.genes) %>% inner_join(av.random) %>% inner_join(av.genome)



sum.genes <- right_join(sum.genes,select(dtg,chr,ATG_ID))

