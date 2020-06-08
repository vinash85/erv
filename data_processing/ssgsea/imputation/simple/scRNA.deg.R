
library(data.table)
library(magrittr)
library(Matrix)
library(parallel)
library(Seurat)
library(avinash.scRNA)
library(tidyverse)
# list1 = list.files(full.names =T, recursive =T, pattern = "*.rds", path = "/liulab/zzeng/Tigger/static/data_computed/scICB/")
list2 = list.files(full.names =T, recursive =T, pattern = "*.rds", path = "/liulab/zzeng/Tigger/static/data_private_computed/scICB/")
scrna.seurat.list = c(list2)
scrna.seurat.meta = lapply(scrna.seurat.list, function(tt) {
	type = gsub(".rds", basename(tt),replacement="")
	database = basename(dirname(tt))
	c(fid=tt, type=type, database=database)
}) %>% 
do.call(rbind,.) %>% as.data.table() %>%
.[,lab:=paste(type, database, sep=":")]
dim(scrna.seurat.meta)
length(unique(scrna.seurat.meta$lab))
# scrna.seurat.meta 

mydeg.tryCatch = function(tt,...){
	tryCatch(mydeg(tt,...),
		error=function(e) NA)

}

scrna.dataset.list = list()
for (ii in seq(nrow(scrna.seurat.meta))) {
	sco = readRDS(file=scrna.seurat.meta[ii]$fid)
	sco$patient.name = sco$sample 
	scrna.dataset.list[[scrna.seurat.meta[ii]$lab]] =  mydeg.tryCatch(sco)
}
scrna.dataset.list = scrna.dataset.list[!is.na(scrna.dataset.list)]
all.genes = scrna.dataset.list %>% sapply(., function(tt) as.character(tt[["gene"]])) %>% unlist %>% unique
scrna.dataset.deg.mat = matrix(NA, nrow=length(all.genes), ncol=length(scrna.dataset.list),
	dimnames =list(all.genes, names(scrna.dataset.list)))
for (ii in seq(length(scrna.dataset.list))) {
	dt.curr = scrna.dataset.list[[ii]]
	scrna.dataset.deg.mat[dt.curr$gene,ii] = dt.curr$stat
} 
scrna.cor.dataset = cor(scrna.dataset.deg.mat, method="spearman", use="pairwise.complete.obs")
which(scrna.cor.dataset < -0.07, arr.ind = T)



## gse13281
sco = readRDS(file="~/liulab_home/data/single_cell/GSE145281/seurat.RDS") 
sco$patient.name = sco$patient
scrna.dataset.list[["GSE120909.all"]] =  mydeg(sco, 
  responder.string=1, 
  nonresponder.string=0)
sco.cd45 = sco[,sco$seurat_clusters%in%c(0,2,4,6,7,10,13,16,12)]
scrna.dataset.list[["GSE120909.cd45+"]] =  mydeg(sco.cd45, 
  responder.string=1, 
  nonresponder.string=0)

sco.cd45 = sco[,sco$seurat_clusters %in% setdiff(1:16, c(15, 0,2,4,6,7,10,13,16,12))]
scrna.dataset.list[["GSE120909.cd45-"]] =  mydeg(sco.cd45, 
  responder.string=1, 
  nonresponder.string=0)


# gse13281.tcell.deseq = export.rdata(".figs/gse145281/deseq.dt.RData")

## include lee data 
sco = readRDS(file="/liulab/asahu/data/single_cell/Lee_data/lee.seurat.cd3.RDS")
scrna.dataset.list[["Shipp_CD3_positive_pre"]] =  mydeg(sco[,(sco$response %in% c("CR", "PR", "PD")) & sco$Treatment.Cycle=="C1D1"])
scrna.dataset.list[["Shipp_CD3_positive_post"]] =  mydeg(sco[,(sco$response %in% c("CR", "PR", "PD")) &  sco$Treatment.Cycle=="C4D1"])

sco = readRDS(file="/liulab/asahu/data/single_cell/Lee_data/lee.seurat.cd3.neg.RDS")
scrna.dataset.list[["Shipp_CD3_negative_pre"]] =  mydeg(sco[,(sco$response %in% c("CR", "PR", "PD"))  & sco$Treatment.Cycle=="C1D1"])
scrna.dataset.list[["Shipp_CD3_negative_post"]] =  mydeg(sco[,(sco$response %in% c("CR", "PR", "PD")) & sco$Treatment.Cycle=="C4D1"])

saveRDS(file="~/liulab_home/data/single_cell/scrna.dataset.deg.RData",scrna.dataset.list) 