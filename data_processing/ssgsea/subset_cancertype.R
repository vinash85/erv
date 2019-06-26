setwd("/Users/avi/projects/deepImmune/data/tcga/neoantigen.v2/attention/tcga.imputed/brca")
aa = fread("../dataset_train.txt")
bb = aa[cancertype=="BRCA"]
write.table(file="dataset_train.txt", x = bb,
		row.names = F, col.names =T,  sep="\t", quote=F )
aa = fread("../dataset_val.txt")
bb = aa[cancertype=="BRCA"]
write.table(file="dataset_val.txt", x = bb,
		row.names = F, col.names =T,  sep="\t", quote=F )
