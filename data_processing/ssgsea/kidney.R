## how to model kidney cancer

# 1.  TCGA -- survival 


tcga_phenotype = fread("/homes6/asahu/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/dataset.txt")
tcga_phenotype = fread("dataset_tcga.txt")
col1 = colnames(tcga_phenotype)[52:61]
setnames(tcga_phenotype,52:61, paste0(col1, "new") )
samples_name = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2//samples_name.txt")$x
samples_name = fread("samples_name_tcga.txt")$x
# tcga_matched = tcga_phenotype[match(rownames(ref.expression), samples_name)]
tcga_matched = tcga_phenotype


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

# comparison of kirc and breast cancer 

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


factor.common = intersect(kirc.cancer$label, miao.dt$label)
miao.dt1 = miao.dt[match(factor.common, label)]
kirc.cancer1 = kirc.cancer[match(factor.common, label)]
miao.dt1[,effect := ifelse(P< 0.05, coef, 0)]
kirc.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
# avi.dt = data.table(x=miao.dt1$coef, y=kirc.cancer1$coef)
cor.test(miao.dt1$coef, kirc.cancer1$coef)
cor.test(miao.dt1$effect, kirc.cancer1$effect)
# library(ggplot2)
# p = ggplot(aes(x=x, y=y), data=avi.dt) + geom_point()
# ggsave(p, filename="~/project/deeplearning/icb/data/RCC_PD1_Miao/tcga_miao.pdf")

miao.dt.output = miao.dt[grep(label, pattern=".output$")]
miao.dt.output$label = gsub(miao.dt.output$label, pattern=".output$", replacement="")

factor.common = intersect(kirc.cancer$label, miao.dt.output$label)
miao.dt1 = miao.dt.output[match(factor.common, label)]
kirc.cancer1 = kirc.cancer[match(factor.common, label)]
miao.dt1[,effect := ifelse(P< 0.05, coef, 0)]
kirc.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
cor.test(miao.dt1$coef, kirc.cancer1$coef)
cor.test(miao.dt1$effect, kirc.cancer1$effect)

avi.dt  = merge(kirc.cancer1, miao.dt1, by ="label")

p = ggplot(avi.dt, aes(x = coef.x, y = coef.y)) +
  geom_point() +
  # scale_color_manual(values = c("red", "grey")) +

  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x="KIRC-tcga survial impact", y="KIRC-MIAO survival impact") +
  geom_text_repel(
    data = subset(avi.dt,  P.x < 0.01 | P.y < 0.01 ),
    aes(label = label),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1))  

 ggsave(p, file="~/project/deeplearning/icb/data/RCC_PD1_Miao/kindney_tcga_miao.pdf")
 cor.test(avi.dt$coef.x, avi.dt$coef.y)



resp.dt$label = names(resp.associations)
resp.dt = resp.dt[order(V1)]


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


factor.common = intersect(bladder.cancer$label, genentech.dt$label)
genentech.dt1 = genentech.dt[match(factor.common, label)]
bladder.cancer1 = bladder.cancer[match(factor.common, label)]
genentech.dt1[,effect := ifelse(P< 0.05, coef, 0)]
bladder.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
# avi.dt = data.table(x=genentech.dt1$coef, y=bladder.cancer1$coef)
cor.test(genentech.dt1$coef, bladder.cancer1$coef)
cor.test(genentech.dt1$effect, bladder.cancer1$effect)
genentech.dt.output = genentech.dt



genentech.dt.output = genentech.dt[grep(label, pattern=".output$")]
genentech.dt.output$label = gsub(genentech.dt.output$label, pattern=".output$", replacement="")

factor.common = intersect(kirc.cancer$label, genentech.dt.output$label)
genentech.dt1 = genentech.dt.output[match(factor.common, label)]
kirc.cancer1 = kirc.cancer[match(factor.common, label)]
genentech.dt1[,effect := ifelse(P< 0.05, coef, 0)]
kirc.cancer1[,effect := ifelse(P< 0.05, coef, 0)]
cor.test(genentech.dt1$coef, kirc.cancer1$coef)
# cor.test(genentech.dt1$effect, kirc.cancer1$effect)



avi.dt  = merge(bladder.cancer1, genentech.dt1, by ="label")
avi.dt[,PC:=ifelse(grepl(label, pattern="PC"),1,0)]

library(ggplot2)
library(ggrepel)

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



## evaluating weights from attention model. 
 aa = fread("/Users/avi/projects/deepImmune/data/tcga/neoantigen.v2/attention/tcga.imputed/drug/tensorboardLog/20190612-072724/last_val_0.csv")
 # 1. create matrix 
 # 2. match variables and samples 
 # 3. Aggregate across the cancer types --median 
 # 4. correlate 
aa = aa[order(sample_name)]
## 
dataset_val = fread("/Users/avi/projects/deepImmune//data/tcga/neoantigen.v2/attention/tcga.imputed/drug/dataset_val.txt")
colnames(dataset_val)[1:131]
head(colnames(dataset_val))



