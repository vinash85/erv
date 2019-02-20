
file.tpm = "~/Dropbox/project/code/deeplearning/icb/risk/genentech.data/data/Genetech_expression_TPM.txt"
clin.file = "~/Dropbox/project/code/deeplearning/icb/risk/genentech.data/data/IMvigor210.clinical"
geneset.file = "~/Dropbox/project/code/deeplearning/icb/risk/genentech.data/geneset.txt"
output.dir = "~/Dropbox/project/code/deeplearning/icb//data/genentech.tpm/"
dir.create(output.dir)
tpm <- fread(file.tpm)
tpm.mat = t(as.matrix(tpm[,-1, with=F]))
colnames(tpm.mat) = tpm$V1
pd_genes <- c("CD274", "PDCD1")
genes.tr <- fread(geneset.file)$x 
genes <- c(genes.tr,pd_genes)
tpm.mat <- tpm.mat[,(colnames(tpm.mat) %in% genes.tr)]

clin <- fread(clin.file)
dataset_phenotype = clin
setnames(dataset_phenotype,1,"patient.name")
patient.name = rownames(tpm.mat) 

common.patients = intersect(patient.name, dataset_phenotype$patient.name)

phenotype_sel = dataset_phenotype[match(common.patients, dataset_phenotype$patient.name)]

    phenotype_mat =  as.matrix(phenotype_sel[, 2:4,with=F])

dataset_sel = tpm.mat[match(common.patients, patient.name), ] 

######
# correcting for response
#########

sort(phenotype_sel[is.na(Response) & (Event == 1)]$OS)
phenotype_sel.mod = phenotype_sel
phenotype_sel.mod[is.na(Response) & (Event == 1) & (OS < 3)]$Response = 0
phenotype_sel.mod[is.na(Response) & (OS > 7)]$Response = 1


################
# hierarchical clustering
################

mydata =  dataset_sel.mod
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")



rand_inx = sample(nrow(dataset_sel))
dataset_sel_shuffle = dataset_sel[rand_inx,hc$order]
phenotype_mat =  as.matrix(phenotype_sel.mod[, 2:4,with=F])
phenotype_mat_shuffle = phenotype_mat[rand_inx,]
train.inx = 1:ceiling(.5 * nrow(dataset_sel))
val.inx = ceiling(.5 * nrow(dataset_sel)):nrow(dataset_sel)

write.table(file=paste0(output.dir, "/ssgsea_train.txt"),x = dataset_sel_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/phenotype_train.txt"),x = phenotype_mat_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )

write.table(file=paste0(output.dir, "/ssgsea_val.txt"),x = dataset_sel_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/phenotype_val.txt"),x = phenotype_mat_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )



################
# PCA
################

# standardise
# minmax <- function(x) (x - min(x))/(max(x) - min(x))
# x_train <- apply(dataset_sel, 2, minmax)
x_train = dataset_sel
pca <- prcomp(t(x_train), center = TRUE, scale. = TRUE)

pca_dataset = pca$rotation

rand_inx = sample(nrow(pca_dataset))
pca_dataset_shuffle = pca_dataset[rand_inx,]
phenotype_mat =  as.matrix(phenotype_sel.mod[, 2:4,with=F])
phenotype_mat_shuffle = phenotype_mat[rand_inx,]
train.inx = 1:ceiling(.5 * nrow(pca_dataset))
val.inx = ceiling(.5 * nrow(pca_dataset)):nrow(pca_dataset)
output.dir = "~/Dropbox/project/code/deeplearning/icb//data/genentech.pca/" 
write.table(file=paste0(output.dir, "/ssgsea_train.txt"),x = pca_dataset_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/phenotype_train.txt"),x = phenotype_mat_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )

write.table(file=paste0(output.dir, "/ssgsea_val.txt"),x = pca_dataset_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/phenotype_val.txt"),x = phenotype_mat_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
