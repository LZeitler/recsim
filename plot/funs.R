guessstart <- function(data){
    dt <- data
    data.frame(gen=rep(-1,max(dt$rep)+1),
               freq=rep(round(mean(filter(dt,gen==0)$freq),1),max(dt$rep)+1),
               rep=c(0:max(dt$rep)+1))
}

dip.pops <- c('BEL','BIH','CRO','FOJ','GOR','HNE','HNI','KZL','MIE','PRE','RZA','SNO','SZI','TRD','VEL','VID')
