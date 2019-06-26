## how to model kidney cancer

# 1.  TCGA -- survival 
library(survival)

setwd("~/Dropbox/project/code/deeplearning/icb/results/jun10/")

source("~/Dropbox/project/code/deeplearning/icb/deepImmune/source.R")
tcga_phenotype = fread("/homes6/asahu/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/dataset.txt")

tcga_phenotype = fread("dataset_tcga.txt")
col1 = colnames(tcga_phenotype)[52:61]
setnames(tcga_phenotype,52:61, paste0(col1, "new") )
samples_name = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2//samples_name.txt")$x
samples_name = fread("samples_name_tcga.txt")$x
# tcga_matched = tcga_phenotype[match(rownames(ref.expression), samples_name)]
tcga_matched = tcga_phenotype

plot.heatmap = function(dat, filename, height = 10, width =7){
  hc = hclust(as.dist(1-cor(dat, method="spearman", use="pairwise.complete.obs")), method="complete")
  hr = hclust(as.dist(1-cor(t(dat), method="spearman", use="pairwise.complete.obs")), method="complete")

  require(heatmap3)
  pdf( filename, height = height, width =width)

  heatmap3(dat, Rowv=as.dendrogram(hr),  Colv=as.dendrogram(hc), scale="none", balanceColor=T, showRowDendro=T ,   showColDendro=T, cexRow = .5, cexCol = 1)

  dev.off()
}


cox.association = function(survival, cov){
	require(survival)
	cov = unlist(cov)
	cov = qnorm.array(cov)
	tryCatch({
		dt1 = data.table(survival, cov=cov)
		setnames(dt1,1:2,c("time", "status"))
		cox.out = coxph(Surv(time,status) ~., data=dt1)
		aa  = summary(cox.out)
		aa$coefficients["cov",]
		},
		error = function(e) NA
		)

	}

response.association = function(response, cov){
	cov = unlist(cov)
	tryCatch(
	{
	aa = wilcox.test(cov[response=="R"], cov[response=="NR"])
	c(aa$p.value, mean(cov[response=="R"], na.rm=T)/( mean(cov[response=="NR"],na.rm=T) + 1E-3))
		},
		error = function(e) NA
		)
		
	}




# comparison of kirc and breast cancer 
cancer.types = c("breast", "bladder", "melanoma", "kidney", "HNSC", "LUSC", "LUAD", "PRAD")
cancer.type.order = t(t(sort(table(tcga_phenotype$cancertype),decreasing=T)))

# cancer.types = cancer.codes = unique(c("BRCA", "BLCA", "SKCM", "KIRC", "HNSC", "LUSC", "LUAD", "PRAD"))
cancer.types = cancer.codes = rownames(cancer.type.order)
names(cancer.codes) = cancer.types
type.combo = combn(cancer.types, 2)

surv.cancers = list()

for (cancer.type in cancer.types) {
tcga.rcc = tcga_phenotype[cancertype==cancer.codes[cancer.type],] 
surv.inx = 63:64
tcga.pheno.mat = tcga.rcc[,-surv.inx,with=F]
tcga.survival = tcga.rcc[,surv.inx,with=F]

surv.associations = mclapply(seq(2,ncol(tcga.pheno.mat)), function(tt) cox.association(tcga.survival, tcga.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(tcga.pheno.mat)[seq(2,ncol(tcga.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
# surv.dt = surv.dt[order(P)]
surv.cancers[[cancer.type]] = surv.dt 
print(cancer.type)

}


## create heatmap 
label.order = surv.cancers[[1]]$label
surv.cancers.coefs = lapply(surv.cancers, function(tt) 
  {
    tt = tt[match(label.order, label)]
    aa = ifelse(tt$coef > 0 , -log10(tt$P), log10(tt$P))
    # aa = ifelse(abs(aa) > 5, 5*sign(aa), aa)
    aa
  }
    )
surv.cancers.coefs = do.call(cbind, surv.cancers.coefs)
rownames(surv.cancers.coefs) = label.order
surv.cancers.dt = surv.cancers.coefs[rowSums(is.na(surv.cancers.coefs)) < 30, ] 
 plot.heatmap(dat=surv.cancers.dt , filename = "~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/pancancer_immune_factors_heatmap_pvalue.pdf", height = 20, width =10)


surv.cancers.coefs = lapply(surv.cancers, function(tt) 
  {
    tt = tt[match(label.order, label)]
    aa = ifelse(tt$P < 0.01, 0, tt$coef)
    aa = ifelse(abs(aa) > 1.5, 1.5*sign(aa), aa)
    aa
  }
    )
surv.cancers.coefs = do.call(cbind, surv.cancers.coefs)
rownames(surv.cancers.coefs) = label.order
surv.cancers.dt = surv.cancers.coefs[rowSums(is.na(surv.cancers.coefs)) < 30, ] 
 plot.heatmap(dat=surv.cancers.dt , filename = "~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/pancancer_immune_factors_heatmap_coef.pdf", height = 20, width =10)



# for (inx in seq(ncol(type.combo))) {
library(parallel)
out =lapply( seq(ncol(type.combo)),  function(inx) {
  combo.curr = type.combo[,inx]

  avi.dt =merge(surv.cancers[[combo.curr[1]]], surv.cancers[[combo.curr[2]]], by='label')
  avi.dt[,coef.x:=ifelse(coef.x > 0.5 , 0.5, ifelse(coef.x < -0.5 , -0.5, coef.x))]
  avi.dt[,coef.y:=ifelse(coef.y > 0.5 , 0.5, ifelse(coef.y < -0.5 , -0.5, coef.y))]

  p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x= sprintf("%s (TCGA) survial impact", combo.curr[1]), 
   y=sprintf("%s (TCGA) survial impact", combo.curr[2])) +
  geom_text_repel(
    data = subset(avi.dt,  (P.x < 0.01 | P.y < 0.01)),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.5, 0.5))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

filename=sprintf("~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/%s_%s_tcga.pdf", combo.curr[1], combo.curr[2])
 ggsave(p, file=filename)

})


## ICB cohorts

# 2. RCC miao et. al.
# icb.phenotype = fread("~/project/deeplearning/icb/data/RCC_PD1_Miao/val_prediction.csv")
# icb.samples = fread("~/project/deeplearning/icb/data/RCC_PD1_Miao/sample_names.txt")
# icb.dataset = fread("~/project/deeplearning/icb/data/RCC_PD1_Miao/dataset.txt")

icb.phenotype = fread("val_prediction_miao.txt")
icb.samples = fread("sample_names_miao.txt")
icb.dataset = fread("dataset_miao.txt")

icb.phenotype = icb.phenotype[unlist(icb.phenotype$sample_name) +1]
icb.response = ifelse(icb.dataset$imputed.response==1 , "R", "NR")
icb.survival = icb.dataset[,.(survive, vital_status)]

icb.pheno.mat = icb.phenotype[,-1,with=F]
resp.associations = apply(icb.pheno.mat, 2, function(tt) response.association(icb.response, tt))
resp.dt = data.table(do.call(rbind, resp.associations))
resp.dt$label = names(resp.associations)
resp.dt = resp.dt[order(V1)]

surv.associations = mclapply(seq(1,ncol(icb.pheno.mat)), function(tt) cox.association(icb.survival, icb.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(icb.pheno.mat)[seq(1,ncol(icb.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
# surv.dt = surv.dt[order(P)]
miao.dt = surv.dt


factor.common = intersect(surv.cancers[["KIRC"]]$label, miao.dt$label)
miao.dt1 = miao.dt[match(factor.common, label)]
kirc.cancer1 = kirc.cancer[match(factor.common, label)]
miao.dt1[,effect := ifelse(P< 0.05, coef, 0)]
kirc.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
cor.test(miao.dt1$coef, kirc.cancer1$coef)
cor.test(miao.dt1$effect, kirc.cancer1$effect)

miao.dt.output = miao.dt[grep(label, pattern=".output$")]
miao.dt.output$label = gsub(miao.dt.output$label, pattern=".output$", replacement="")

factor.common = intersect(surv.cancers[["KIRC"]]$label, miao.dt.output$label)
miao.dt1 = miao.dt.output[match(factor.common, label)]
# kirc.cancer1 = kirc.cancer[match(factor.common, label)]
# miao.dt1[,effect := ifelse(P< 0.05, coef, 0)]
# kirc.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
# cor.test(miao.dt1$coef, kirc.cancer1$coef)
# cor.test(miao.dt1$effect, kirc.cancer1$effect)
avi.dt  = merge(surv.cancers[["KIRC"]], miao.dt1, by ="label")

  avi.dt[,coef.x:=ifelse(coef.x > 1 , 1, ifelse(coef.x < -1 , -1, coef.x))]
  avi.dt[,coef.y:=ifelse(coef.y > 1 , 1, ifelse(coef.y < -1 , -1, coef.y))]
library(ggrepel)
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="KIRC-tcga survial impact", y="MIAO survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  P.x < 0.01 | P.y < 0.01 ),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.5, 0.5))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


dir.create("~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/icb_cohort_comparsion")
 ggsave(p, file="~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/icb_cohort_comparsion/KIRC_miao.pdf")
 # cor.test(avi.dt$coef.x, avi.dt$coef.y)



## genentech data 

## file transfered from  local mac ../data/tcga/neoantigen.v2/attention/genentech.imputed/val_prediction.csv
icb.phenotype = fread("~/genentech_val.csv")
icb.phenotype = fread("genentech_val.csv")
icb.input = fread("/Users/avi/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/genentech.imputed/dataset.txt")
icb.survival = icb.phenotype[,.(survive, vital_status)]
surv.associations = mclapply(seq(1,ncol(icb.pheno.mat)), function(tt) cox.association(icb.survival, icb.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(icb.pheno.mat)[seq(1,ncol(icb.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
# surv.dt = surv.dt[order(P)]
genentech.dt = surv.dt


factor.common = intersect(surv.cancers[["BLCA"]]$label, genentech.dt$label)
genentech.dt1 = genentech.dt[match(factor.common, label)]
# surv.cancers[["BLCA"]] = surv.cancers[["BLCA"]][match(factor.common, label)]
# genentech.dt1[,effect := ifelse(P< 0.05, coef, 0)]
# surv.cancers[["BLCA"]][,effect := ifelse(P< 0.05, coef, 0)]
# avi.dt = data.table(x=genentech.dt1$coef, y=surv.cancers[["BLCA"]]1$coef)
cor.test(genentech.dt1$coef, surv.cancers[["BLCA"]]$coef)
cor.test(genentech.dt1$effect, surv.cancers[["BLCA"]]$effect)
genentech.dt.output = genentech.dt



genentech.dt.output = genentech.dt[grep(label, pattern=".output$")]
genentech.dt.output$label = gsub(genentech.dt.output$label, pattern=".output$", replacement="")

# factor.common = intersect(kirc.cancer$label, genentech.dt.output$label)
# genentech.dt1 = genentech.dt.output[match(factor.common, label)]
# kirc.cancer1 = kirc.cancer[match(factor.common, label)]
# genentech.dt1[,effect := ifelse(P< 0.05, coef, 0)]
# kirc.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
# cor.test(genentech.dt1$coef, kirc.cancer1$coef)
# cor.test(genentech.dt1$effect, kirc.cancer1$effect)



avi.dt  = merge(surv.cancers[["BLCA"]], genentech.dt.output, by ="label")

  avi.dt[,coef.x:=ifelse(coef.x > 1 , 1, ifelse(coef.x < -1 , -1, coef.x))]
  avi.dt[,coef.y:=ifelse(coef.y > 1 , 1, ifelse(coef.y < -1 , -1, coef.y))]
library(ggrepel)
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="Bladder-tcga survival impact", y="Bladder-Genentech survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  P.x < 0.1 | P.y < 0.1 ),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.5, 0.5))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


 ggsave(p, file="~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/icb_cohort_comparsion/BLCA_genentech.pdf")



############
### old 
############

tcga.rcc = tcga_phenotype[cancertype=="BRCA",] 
surv.inx = 63:64
tcga.pheno.mat = tcga.rcc[,-surv.inx,with=F]
tcga.survival = tcga.rcc[,surv.inx,with=F]

library(parallel)
library(survival)
surv.associations = mclapply(seq(2,ncol(tcga.pheno.mat)), function(tt) cox.association(tcga.survival, tcga.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(tcga.pheno.mat)[seq(2,ncol(tcga.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
surv.dt = surv.dt[order(P)]

breast.cancer = surv.dt 

## multivariate analysis 
cols = paste0("PC", c(11,15,48))
temp.dt = tcga.pheno.mat[ , (cols) := lapply(.SD, qnorm.array), .SDcols = cols]
comb.dt = cbind( temp.dt, tcga.survival)
suma = coxph(Surv(survive, vital_status) ~ PC11 + PC15 + PC48, comb.dt)
comb.dt[,surv.pred:=.2*PC15 -.3*PC48]
suma = coxph(Surv(survive, vital_status) ~ surv.pred, comb.dt)
concordance(suma)


###########
##bcla
#######
tcga.rcc = tcga_phenotype[cancertype=="BLCA",] 
surv.inx = 63:64
tcga.pheno.mat = tcga.rcc[,-surv.inx,with=F]
tcga.survival = tcga.rcc[,surv.inx,with=F]

surv.associations = mclapply(seq(2,ncol(tcga.pheno.mat)), function(tt) cox.association(tcga.survival, tcga.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(tcga.pheno.mat)[seq(2,ncol(tcga.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
surv.dt = surv.dt[order(P)]

bladder.cancer = surv.dt 


tcga.rcc = tcga_phenotype[cancertype=="SKCM",] 
surv.inx = 63:64
tcga.pheno.mat = tcga.rcc[,-surv.inx,with=F]
tcga.survival = tcga.rcc[,surv.inx,with=F]

surv.associations = mclapply(seq(2,ncol(tcga.pheno.mat)), function(tt) cox.association(tcga.survival, tcga.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(tcga.pheno.mat)[seq(2,ncol(tcga.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
surv.dt = surv.dt[order(P)]

melanoma.cancer = surv.dt 


tcga.rcc = tcga_phenotype[cancertype=="KIRC",] 
surv.inx = 63:64
tcga.pheno.mat = tcga.rcc[,-surv.inx,with=F]
tcga.survival = tcga.rcc[,surv.inx,with=F]

surv.associations = mclapply(seq(2,ncol(tcga.pheno.mat)), function(tt) cox.association(tcga.survival, tcga.pheno.mat[,tt,with=F]), mc.cores=32)

surv.dt = do.call(rbind, surv.associations)
surv.dt = data.table(surv.dt)
surv.dt$label = colnames(tcga.pheno.mat)[seq(2,ncol(tcga.pheno.mat))]
setnames(surv.dt, 2:5, c("exp", "se", "z", "P"))
surv.dt = surv.dt[order(P)]
kirc.cancer = surv.dt
kidney.cancer = surv.dt 



avi.dt  = merge(miao.dt, genentech.dt.output, by ="label")

  avi.dt[,coef.x:=ifelse(coef.x > 1 , 1, ifelse(coef.x < -1 , -1, coef.x))]
  avi.dt[,coef.y:=ifelse(coef.y > 1 , 1, ifelse(coef.y < -1 , -1, coef.y))]
library(ggrepel)
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="Miao survival impact", y="Bladder-Genentech survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  P.x < 0.1 | P.y < 0.1 ),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.5, 0.5))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


 ggsave(p, file="~/Dropbox/project/code/deeplearning/icb/results/jun10/survival_comparsion/icb_cohort_comparsion/miao_genentech.pdf")






p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  # geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="Bladder-tcga survival impact", y="Bladder-Genentech survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  (P.x < 0.1 | P.y < 0.1) & PC ==0 ),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.25, 0.25)) + 
  # ylim(c(-0.25, 0.25))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)



 ggsave(p, file="~/project/deeplearning/icb/data/RCC_PD1_Miao/bladder_tcga_genentech.pdf")
 ggsave(p, file="bladder_tcga_genentech.pdf")
 cor.test(avi.dt$coef.x, avi.dt$coef.y)


## genentech and miao et. al. 
avi.dt  = merge(miao.dt1, genentech.dt1, by ="label")
avi.dt[,PC:=ifelse(grepl(label, pattern="PC"),1,0)]
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="KIRC-MIAO survival impact", y="Bladder-Genentech survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  (P.x < 0.06 | P.y < 0.06) & PC ==0 ),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.25, 0.25)) + 
  # ylim(c(-0.25, 0.25))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)



 ggsave(p, file="~/project/deeplearning/icb/data/RCC_PD1_Miao/miao_genentech.pdf")
 cor.test(avi.dt$coef.x, avi.dt$coef.y)







 # cor.test(avi.dt$coef.x, avi.dt$coef.y)





avi.dt =merge(kirc.cancer, breast.cancer, by='label')
library(ggrepel)
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  geom_point() +
  # scale_color_manual(values = c("red", "grey")) +

  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="KIRC survial impact", y="BRCA survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  P.x < 0.01 | P.y < 0.01 ),
    aes(label = label),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  xlim(c(-0.5, 0.5)) + 
  ylim(c(-0.5, 0.5))  

 ggsave(p, file="~/project/deeplearning/icb/data/RCC_PD1_Miao/kindney_breast_tcga.pdf")
 cor.test(avi.dt$coef.x, avi.dt$coef.y)



avi.dt =merge(bladder.cancer, melanoma.cancer, by='label')
avi.dt = avi.dt[!(is.na(coef.x)|is.na(coef.y))]
avi.dt[,PC:=ifelse(grepl(label, pattern="PC"),1,0)]
library(ggrepel)
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="BLCA-tcga survival impact", y="SKCM-tcga survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  (P.x < 0.01 | P.y < 0.01) & PC ==0 ),
    aes(x = coef.x, y = coef.y, label = label),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  xlim(c(-0.5, 0.5)) + 
  ylim(c(-0.5, 0.5))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)



 ggsave(p, file="~/project/deeplearning/icb/data/RCC_PD1_Miao/Bladder_melanoma_tcga.pdf")
 cor.test(avi.dt$coef.x, avi.dt$coef.y, method="spearman")


# ggplot(simdata) +
#   theme_geometry(simdata$x, simdata$y) +
#   geom_point(aes(x = x, y = y), size = 3, color = "red") + 
#   ggtitle("More geometric example")



avi.dt =merge(kirc.cancer, melanoma.cancer, by='label')
library(ggrepel)
p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  geom_point() +
  # scale_color_manual(values = c("red", "grey")) +

  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="KIRC-tcga survival impact", y="SKCM-tcga survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  P.x < 0.01 | P.y < 0.01 ),
    aes(label = label),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  xlim(c(-0.5, 0.5)) + 
  ylim(c(-0.5, 0.5))  

 ggsave(p, file="~/project/deeplearning/icb/data/RCC_PD1_Miao/kindney_melanoma_tcga.pdf")
 cor.test(avi.dt$coef.x, avi.dt$coef.y)





## evaluating weights from attention model. 
 # aa = fread("/Users/avi/projects/deepImmune/data/tcga/neoantigen.v2/attention/tcga.imputed/drug/tensorboardLog/20190612-072724/last_val_0.csv")
 aa = fread("~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/brca/tensorboardLog/20190622-083937/best_val_0.csv")
 # 1. create matrix 
 # 2. match variables and samples 
 # 3. Aggregate across the cancer types --median 
 # 4. correlate 
aa = aa[order(sample_name)]
## 
# dataset_val = fread("/Users/avi/projects/deepImmune//data/tcga/neoantigen.v2/attention/tcga.imputed/drug/dataset_val.txt")
dataset_val = fread("~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/brca/dataset_val.txt")
aa.inx = unlist(c(list(seq(2,62)), list(seq(66, 137)) , list(seq(138,169))))
input.names = colnames(dataset_val)[aa.inx][1:131]

aa.mat = aa[, seq(267-130, 267), with=F]
setnames(aa.mat, 1:131, input.names)
aa.mat$cancertype = dataset_val$cancertype

melanoma.attn = aa.mat[cancertype=="SKCM"]
attn.mean = apply(melanoma.attn[,1:131,with=F], 2,median)
melanoma.cancer$attn.mean = attn.mean[melanoma.cancer$label]
cor.test(melanoma.cancer[P<0.01]$coef, melanoma.cancer[P<0.01]$attn.mean, method="spearman")


p = ggplot(melanoma.cancer, aes(x = coef, y = attn.mean)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  # geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="SKCM-tcga survival impact", y="SKCM mean attention") +
  geom_text_repel(
    data = melanoma.cancer,
    aes(x = coef, y = attn.mean, label = label),
    size = 1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.25, 0.25)) + 
  # ylim(c(-0.25, 0.25))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

 ggsave(p, file="skcm_tcga_attention.pdf")


## breast cancer



breast.attn = aa.mat[cancertype=="BRCA"]
attn.mean = apply(breast.attn[,1:131,with=F], 2,median)
breast.cancer$attn.mean = attn.mean[breast.cancer$label]
cor.test(breast.cancer[P<0.01]$coef, breast.cancer[P<0.01]$attn.mean, method="spearman")


p = ggplot(breast.cancer, aes(x = coef, y = attn.mean)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  # geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="BRCA-tcga survival impact", y="BRCA mean attention") +
  geom_text_repel(
    data = breast.cancer,
    aes(x = coef, y = attn.mean, label = label),
    size = 1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.25, 0.25))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

 ggsave(p, file="brca_tcga_attnetion.pdf")



## comparison of attnetion between brca and melanoma

breast.attn.mean = apply(breast.attn[,1:131,with=F], 2,median)
melanoma.attn.mean = apply(melanoma.attn[,1:131,with=F], 2,median)

cor.test(breast.attn.mean, melanoma.attn.mean) 
avi.dt = data.table(breast.attn = breast.attn.mean, melanoma.attn = melanoma.attn.mean, label = names(breast.attn.mean))

p = ggplot(avi.dt, aes(x = breast.attn, y = melanoma.attn)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  # geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="BRCA mean attention", y="Melanoma mean attention") +
  geom_text_repel(
    data = avi.dt,
    aes(x = breast.attn, y = melanoma.attn, label = label),
    size = 1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  # xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.25, 0.25))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

 ggsave(p, file="brca_melanoma_attnetion.pdf")



## 3 variables debug specifically in breast cancer
 # aa = fread("~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/brca/tensorboardLog/three_var_20190622-191219/last_val_0.csv")
 aa = fread("~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/brca//val_prediction.csv")

suma = coxph(Surv(survive, vital_status) ~ survive.output, aa)
concordance(suma)
aa = aa[order(sample_name)]
dataset_val = fread("~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/brca/dataset_val.txt")
# dataset_val = fread("~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/dataset.txt")
dataset_val = dataset_val[cancertype=="BRCA"]
aa.inx = unlist(c(list(seq(2,62)), list(seq(66, 137)) , list(seq(138,169))))
input.names = colnames(dataset_val)[aa.inx][c(11,15, 48)]

aa.mat = aa[, seq( ncol(aa) - length(input.names) + 1, ncol(aa)), with=F]
setnames(aa.mat, input.names)
aa.mat$cancertype = dataset_val$cancertype

breast.attn = aa.mat[cancertype=="BRCA"]
attn.mean = apply(breast.attn[,-ncol(aa.mat),with=F], 2,median)
breast.cancer$attn.mean = attn.mean[breast.cancer$label]

## check survival make sense
temp1 = dataset_val
temp1[,(cols):=lapply(.SD, qnorm.array),  .SDcols = cols]
breast.prod = breast.attn[, .(PC11, PC15, PC48)] *  dataset_val[, .(PC11, PC15, PC48)]
score = rowSums(breast.prod)

library(survival)
dt = data.table(survive=dataset_val$survive, vital_status=dataset_val$vital_status,qnorm.array(score))
suma = coxph(Surv(survive, vital_status) ~ score, dt)
concordance(suma)

suma = coxph(Surv(dt$survive, dt$vital_status) ~ aa$survive.output, dt)
concordance(suma)
# dt1 = dataset_val[,.(survive, vital_status, PC11, PC15, PC48)], breast.attn[, .(PC11, PC15, PC48)])

suma = coxph(Surv(survive, vital_status) ~ PC11 + PC15 + PC48, dataset_val)



cor.test(breast.cancer[P<0.01]$coef, breast.cancer[P<0.01]$attn.mean, method="spearman")


p = ggplot(breast.cancer, aes(x = coef, y = attn.mean)) +
  # theme_geometry(avi.dt$coef.x, avi.dt$coef.y) +
  geom_point() +
  # geom_smooth(method = lm, se = FALSE) + 
  # scale_color_manual(values = c("red", "grey")) +

  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="BRCA-tcga survival impact", y="BRCA mean attention") +
  geom_text_repel(
    data = breast.cancer,
    aes(x = coef, y = attn.mean, label = label),
    size = 1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  xlim(c(-0.5, 0.5)) + 
  # ylim(c(-0.25, 0.25))  + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

 ggsave(p, file="brca_tcga_attnetion.pdf")




