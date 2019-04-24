load("temp.RData")
hist(aa, 100, col=rgb(1,0,0,0.5), xlim=c(0,27))
hist(bb, 100, col=rgb(0,0,1,0.5), add=T)

aa = c(rep(1,100), rep(0,100))
dens.obs = density(aa, adjust=0.1)
resample.obs = sample(dens.obs$x, 1000, replace=TRUE, prob=dens.obs$y)
dat = data.frame(value=c(aa,resample.obs), 
                 group=rep(c("Observed","Modeled"), c(length(aa),length(resample.obs))))

library(ggplot2)
ggplot(dat, aes(value, fill=group, colour=group)) +
  stat_ecdf(geom="step") +
  theme_bw()


xx1 = match.distribution.zeros( c(rep(0,100), rep(1,100)), aa)
xx1 = match.distribution.zeros(aa, c(rep(0,100)))

aa = match.expression.distribution(dataset_ssgsea_sel[,1:100], ref.expression.cancertype[,1:100])
aa = match.distribution.zeros( dataset_ssgsea_sel[,"DPM1"], ref.expression.cancertype[,"DPM1"])

# association of pi3k pathway and NK cell activity 
# cor.test(dataset_ssgsea_sel.back[,"PIK3CA"], dataset$MLH1)

