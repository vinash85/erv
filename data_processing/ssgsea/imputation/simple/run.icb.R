
## convert to Seurat object and preprocess
library(data.table)
library(magrittr)
library(ggplot2)
library(Matrix)
library(Seurat)
library(cowplot)
library(uwot)
library(parallel)
library(EnhancedVolcano)
library(tidyverse)
options(error=recover)
library(matrixTests)



expression.files = list.files(path = "~/liulab_home/data/ssgsea/xiaoman/icb/expression/", full.names = T)

create.deepImmune.input <- function(output.dir, expression.file,  genes.curr) {
  aa = fread(expression.file)
  genes.inx = match(genes.curr, aa$Symbol)
  exp.mat = as.matrix(t(aa[genes.inx,-1,with=F]) )
  colnames(exp.mat) = genes.curr
  exp.mat
}
load("~/project/deeplearning/icb/data/tcga/tf/genes.analyzed.RData")
dataset.exp.only =  sapply( expression.files, function(expression.file) create.deepImmune.input(expression.file = expression.file,  genes.curr=genes.analyzed), USE.NAMES = T)
dataset.name = sapply(expression.files, basename)
expression.dataset.patinet.name = lapply(dataset.exp.only, rownames)
names(expression.dataset.patinet.name) = dataset.name
all.expression = do.call(rbind, dataset.exp.only)
dataset.label = sapply(seq(length(dataset.exp.only)), function(tt) rep(tt, nrow(dataset.exp.only[[tt]])))


expression.mat.list  = sapply(expression.files, function(expression.file) {
  aa = fread(expression.file)
  genes.curr= aa$Symbol
  exp.mat = as.matrix(t(aa[,-1,with=F]) )
  colnames(exp.mat) = genes.curr
  exp.mat = apply(exp.mat,2,avinash::znorm)
  exp.mat
}
)
all.genes  = unique(unlist(sapply(expression.mat.list, colnames)))
dataset.name  = sapply(expression.files, basename)
names(expression.mat.list) = dataset.name
expression.merge = do.call(rbind, sapply(expression.mat.list, function(tt) tt[,match(all.genes, colnames(tt))], simplify = F))
nan.inx = apply(expression.merge,2,function(tt) sum(!is.na(tt)))
expression.merge = expression.merge[,nan.inx > 200]
all.genes.sel = all.genes[nan.inx > 200]
colnames(expression.merge) = all.genes.sel



followups = sapply(list.files(path = "~/liulab_home/data/ssgsea/xiaoman/icb/follow_up/", full.names = T), fread, USE.NAMES = TRUE)
dataset.name = sapply(list.files(path = "~/liulab_home/data/ssgsea/xiaoman/icb/expression/", full.names = T), basename)
meta.info.dt = do.call(rbind, lapply(dataset.name, function(tt) data.table(Patient=expression.dataset.patinet.name[[tt]], names=tt)))
meta.info.dt$inx = seq(nrow(meta.info.dt))
meta.info.dt[,label:=paste(Patient, names)]
names(followups) = sapply(names(followups), basename)
# sapply(followups, function(tt) setnames(tt, 1, "Patient"))

icb.response.dt =  do.call(rbind, sapply(names(followups), function(tt) {
  aa = followups[[tt]] 
  setnames(aa,1, "Patient")
  aa[,list(Patient, Response, dataset=tt)]
}, simplify = F))

# pfs.dt =  do.call(rbind, sapply(names(followups), function(tt) {
#   tryCatch(followups[[tt]][,list(Patient=V1,  PFS, PFS.Event, dataset=tt)],
#    error = function(e) NULL)
# }, simplify = F))

# os.dt =  do.call(rbind, sapply(names(followups), function(tt) {
#   tryCatch(followups[[tt]][,list(Patient=V1,  OS, OS.Event, dataset=tt)],
#    error = function(e) NULL)
# }, simplify = F))



icb.response.dt = icb.response.dt[!is.na(Response)]
icb.response.dt[,Response1:=as.integer(ifelse(Response> 0, 1,0))]

expression.merge.match = expression.merge[match(icb.response.dt[,paste(Patient, dataset), ],meta.info.dt$label),]
extract_nmle_table <- function (m1){
  mod = summary(m1)
  beta <- m1$coefficients$fixed[2] #$fixed is not needed
  se <- m1$varFix[2]
  t <- beta/se
  p<- anova(m1)[2,4]
  table=data.frame(cbind(beta,se,t,p))
  return(table)
}
library(nlme)
eval.nlme = function(gene.expression, data.dt){
    # data.dt = data.table with colums  Response1 and dataset (factor to control)
  tryCatch(
  {
    data.dt$col = gene.expression
    data.dt = data.dt[!is.na(col)]
    data.dt = data.dt[dataset %in% names(which(table(data.dt$dataset) > 10))]  # only choose those with at least 10 obervations 
    m1 <- lme(col~ Response1, random=~1|dataset, data=data.dt)
    extract_nmle_table(m1)
  },error=  function(e) rep(NA,4))
}
icb.response.genes = mclapply(seq(ncol(expression.merge.match)), function(tt) eval.nlme(expression.merge.match[,tt], data.dt=icb.response.dt), mc.cores=45)
icb.response.genes.dt =  data.table(do.call(rbind, icb.response.genes))[,genes:=all.genes.sel]
setnames(icb.response.genes.dt, 1:4, c("estimate","se", "z", "P"))
library(magrittr)
icb.response.genes.dt = icb.response.genes.dt %>%
.[!is.na(estimate)] %>%
.[order(P)]  %>% 
.[,effect:=sign(estimate)*ifelse(abs(estimate)-se < 0, 0, abs(estimate)-se)]


expression.merge.match = expression.merge[match(icb.response.dt[,paste(Patient, dataset), ],meta.info.dt$label),]


differential.expression.bulk <- function(expression.mat, icb.response.dt) {

  dataset.inx = unique(icb.response.dt$dataset)
  icb.genes = rownames(expression.mat)
  all.stat = list()
  for (grp in dataset.inx) {
   icb.response.curr = icb.response.dt[(dataset == grp)]
   exp.curr = expression.mat[,icb.response.curr$inx]
   responder.grp = which(icb.response.curr$response==1) 
   nonresponder.grp = which(icb.response.curr$response==0) 
   x = exp.curr[,responder.grp]
   y = exp.curr[, nonresponder.grp]
   aa = matrixTests::row_wilcoxon_twosample(x,y)
   all.stat[[grp]] = data.table(P=aa$pvalue, responder.mean = rowMeans(x,na.rm = T), nonresponder.mean = rowMeans(y,na.rm = T)) %>%
   .[,Padj:=p.adjust(P)] %>%
   .[,type:=grp] %>%
   .[,gene:=icb.genes] 

 }

 all.stat.dt = data.table(do.call(rbind, all.stat))
 all.stat.dt
}

icb.response.dt[,inx:=seq(.N)][,response:=Response1]
bulk.icb.dt = differential.expression.bulk(t(expression.merge.match), icb.response.dt)
bulk.icb.dt = bulk.icb.dt[!is.na(P)]
bulk.combined.dt = rbind(
  bulk.icb.dt %>%
  .[,up_or_down:=ifelse(responder.mean > nonresponder.mean, "upregulated", "downregulated")] %>%
  .[,.(P, up_or_down, Padj, type, gene=gene)],
  icb.response.genes.dt %>%
  .[,type:="combined"] %>%
  .[,up_or_down:=ifelse(estimate > 0, "upregulated", "downregulated")] %>%
  .[,Padj:=p.adjust(P)] %>%
  .[,.(P, up_or_down, Padj, type, gene=genes)])

bulk.combined.dt[,deg.effect:=ifelse(up_or_down=="upregulated",-1, 1) * log(P)]

saveRDS(file="/homes6/asahu/liulab_home/data/immunotherapy-trials/all.bulk.rnaseq.deg.Rds", bulk.combined.dt)
