source('~/code/r/source_me.R')
source('~/pro/recsim/plot/funs.R')
setwd('~/pro/recsim/300_analyses/pca/')

## BiocManager::install("SNPRelate")
library(SNPRelate)
library(ggrepel)



pcaplot <- function(data,ev.x,ev.y,...){
    p <- ggplot(data)
    p <- p + geom_point(size=2,alpha=.7,aes(...)) + labs(x = paste('EV', ev.x,
                                                 round(pca.co$varprop[ev.x]*100,2),"%"),
                                       y = paste('EV', ev.y,
                                                 round(pca.co$varprop[ev.y]*100,2),"%"))
    p <- p + scale_alpha_discrete(range = c(.4, .9))+
        scale_shape_manual(values = c(2,19,21))
    ## p <- p + scale_color_viridis(discrete = T)
    return(p)
}



gdsfile <- '~/pro/recsim/data/DM_HM_BP_BI_merged_fi_dip.gds'
snpgdsVCF2GDS('/nfs/nas22/fs2201/biol_impb_group_evo_euler/data/300/PASS/DM_HM_BP_BI_merged_fi_dip.vcf',
              gdsfile,method='biallelic.only')
snpgdsVCF2GDS('/mnt/evo_euler/data/300/filter/DM_HM_BP_BI_merged_fi_dip.vcf',
              gdsfile,method='biallelic.only') # local

gd <- snpgdsOpen(gdsfile)

quick <- fread('../annotations/300_quicklist.csv',data.table=F)

makepca <- function(gdsfile,...){
    pca.co <- snpgdsPCA(gdsfile, autosome.only=F, maf=.01, missing.rate = .1,...)
    pca.co.tab <- data.frame(sample.id = pca.co$sample.id,
                             pop = as.factor(substr(pca.co$sample.id,1,3)),
                             ev1 = pca.co$eigenvect[,1],
                             ev2 = pca.co$eigenvect[,2],
                             ev3 = pca.co$eigenvect[,3],
                             ev4 = pca.co$eigenvect[,4],
                             ev5 = pca.co$eigenvect[,5],
                             stringsAsFactors = F)
}

pca.co.tab <- makepca(gd)

## fwrite(pca.co.tab,paste0('pca_merged_',format(Sys.time(), "%d%m%y_%H%M"),'.txt'))
## pca.co.tab <- fread('pca_merged_100519_1215.txt',data.table=F)

pca.co.tab <- left_join(pca.co.tab,quick,by=c('pop'='population'))
coords <- pca.co.tab %>% group_by(pop) %>% summarise(m_ev1=mean(ev1),m_ev2=mean(ev2))
pca.co.tab <- left_join(pca.co.tab,coords)
    
(p <- pcaplot(pca.co.tab,1,2,x=ev1,y=ev2,color=class,shape=as.character(ploidy))+
    geom_text_repel(data = coords, aes(m_ev1,m_ev2,label=pop)))

saveplot(p,'../plots/pca_snprelate_dip')

pcadat <- filter(pca.co.tab,!class%in%c('Dinaric_2x','Pannonian_2x','Croatica_2x'))
(p1 <- pcaplot(pcadat,1,2,x=ev1,y=ev2,color=class,shape=as.character(ploidy))+
    geom_text_repel(data = filter(coords,pop%in%pcadat$pop), aes(m_ev1,m_ev2,label=pop)))

saveplot(p1,'../plots/pca_snprelate_dip_core')

## pcadat <- filter(pcadat,ploidy==2)
## pcaplot(pcadat,1,2,x=ev1,y=ev2,color=class)+
##     geom_text_repel(data = filter(coords, coords$pop %in% pcadat$pop),
##                     aes(m_ev1,m_ev2,label=pop))

saveplot(p,'../plots/pca_snprelate')

## PCA  with LD pruning

gdp <- snpgdsLDpruning(gd,ld.threshold = .2, autosome.only = F)
pruned <- unlist(unname(gdp))

pca.pr <- makepca(gd,snp.id=pruned)

pca.pr <- left_join(pca.pr,quick,by=c('pop'='population'))
coords <- pca.pr %>% group_by(pop) %>% summarise(m_ev1=mean(ev1),m_ev2=mean(ev2))
pca.pr <- left_join(pca.pr,coords)

(pp <- pcaplot(pca.pr,1,2,x=ev1,y=ev2,color=class,shape=as.character(ploidy))+
    geom_text_repel(data = coords, aes(m_ev1,m_ev2,label=pop)))

saveplot(pp,'../plots/pca_snprelate_dip_pruned')

pcadat <- filter(pca.pr,!class%in%c('Dinaric_2x','Pannonian_2x','Croatica_2x'))
(pp1 <- pcaplot(pcadat,1,2,x=ev1,y=ev2,color=class,shape=as.character(ploidy))+
    geom_text_repel(data = filter(coords,pop%in%pcadat$pop), aes(m_ev1,m_ev2,label=pop)))

saveplot(pp1,'../plots/pca_snprelate_dip_pruned_core')

## prune before PCA calc
pca.pr <- makepca(gd,snp.id=pruned,sample.id=pcadat$sample.id)

pca.pr <- left_join(pca.pr,quick,by=c('pop'='population'))
coords <- pca.pr %>% group_by(pop) %>% summarise(m_ev1=mean(ev1),m_ev2=mean(ev2))
pca.pr <- left_join(pca.pr,coords)

(pp2 <- pcaplot(pca.pr,1,2,x=ev1,y=ev2,color=class,shape=as.character(ploidy))+
    geom_text_repel(data = coords, aes(m_ev1,m_ev2,label=pop)))

saveplot(pp2,'../plots/pca_snprelate_dip_pruned_after_core')

