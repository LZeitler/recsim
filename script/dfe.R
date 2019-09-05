library(dplyr)
library(ggplot2)
library(cowplot)

plotGamma <- function(shape=2, rate=0.5, to=0.99, p=c(0.1, 0.9), cex=1, ...){
  to <- qgamma(p=to, shape=shape, rate=rate)
  curve(dgamma(x, shape, rate), from=0, to=to, n=500, type="l", 
        main=sprintf("gamma(x, shape=%1.2f, rate=%1.2f)", shape, rate),
        bty="n", xaxs="i", yaxs="i", col="blue", xlab="", ylab="", 
        las=1, lwd=2, cex=cex, cex.axis=cex, cex.main=cex, ...)
  gx <- qgamma(p=p,  shape=shape, rate=rate)
  gy <- dgamma(x=gx, shape=shape, rate=rate)
  for(i in seq_along(p)) { lines(x=rep(gx[i], 2), y=c(0, gy[i]), col="blue") }
  for(i in seq_along(p)) { text(x=gx[i], 0, p[i], adj=c(1.1, -0.2), cex=cex) }
}

plotGamma(.3,1)


qgamma(p=.99, shape=.5, rate=1)

dt <- data.frame(fitness=0:100/100,
                 frequency=dgamma(0:100/100, shape = .5))
qplot(dt$fitness,dt$frequency)

myfun <- 'rgamma(10000,shape=.3,rate=.3/.01)*-1  # mean=-0.01'
x <- eval(parse(text=myfun))
## x <- x/max(x)
qplot(x,bins=20)+labs(x="Fitness",y='Frequency',title=myfun)+
    geom_vline(aes(xintercept=mean(x)),color='blue')


myfun <- 'rnorm(1:10000,0,.2)'
x <- eval(parse(text=myfun))
## x <- x/max(x)
qplot(x,bins=20)+labs(x="Fitness",y='Frequency',title=myfun)+
    geom_vline(aes(xintercept=mean(x)),color='blue')


myfun <- 'rexp(10000,1/0.1)  # mean=0.01'           # in eidos: rexp(10000,0.02), 0.02 is the mean, 1/0.02 is
                                        # the rate
x <- eval(parse(text=myfun))
## x <- x/max(x)
qplot(x,bins=20)+labs(x="Fitness",y='Frequency',title=myfun)+
    geom_vline(aes(xintercept=mean(x)),color='blue')

