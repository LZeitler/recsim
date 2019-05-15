source('~/code/r/source_me.R')
source('funs.R')

system('scp euler:pro/recsim/300_analyses/pca/*sc2.eig* ~/pro/recsim/300_analyses/pca/')

val <- fread('../300_analyses/pca/pca_sc2.eigenval',data.table=F)
vec <- fread('../300_analyses/pca/pca_sc2.eigenvec',data.table=F) %>%
    mutate(pop=substr(V1,1,3))

p <- ggplot(vec) +
    geom_point(aes(V3,V4,color=pop))
p

p1 <- ggplot(filter(vec,pop%in%c('KZL','SZI','SNO','RZA'))) +
    geom_point(aes(V3,V4,color=pop))
p1

plot_grid(p,p1)
