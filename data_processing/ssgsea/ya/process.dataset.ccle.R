
# ## Processing of ssgsea data
# Run th file in BCB cluster
# 
# 1. Shape the output ssgsea mat
# 2. Shape the output phenotype
# 3. reorder 
# 4. normalize survival per cancer type
# 5. create 5 NA for continuous output and 5 NA for binary output
# 
# TODO : Write file such sample name is first row and add support to the code

# In[ ]:

#Checking output of pca. prcomp function returns standard deviation (sdev), rotation and loadings.
dataset_ssgsea ="/liulab/asahu/data/ssgsea/xiaoman/ya/gene_13999_expression_crispr_essential_554_cellline.txt"
output.dir = "~/project/deeplearning/icb/data/ya/"

library(data.table)
dir.create(output.dir, showWarnings = FALSE)
cancertype.dt = data.table(ID = dataset_ssgsea$V1, cancertype = as.factor( dataset_ssgsea$cancer_type))
cancertype.pheno = mltools::one_hot(cancertype.dt, sparsifyNAs =T)

dataset = cbind(dataset_ssgsea[,-1,with=F], cancertype.pheno[,-1,with=F])


write.table(file=paste0(output.dir, "/dataset.txt"),x = dataset,
    row.names = F, col.names =T,  sep="\t", quote=F )
# write.table(file=paste0(output.dir, "/dataset_phenotype.txt"),x = phenotype.ext.mat,
    # row.names = F, col.names =T,  sep="\t", quote=F )
rand_inx = sample(nrow(dataset))
dataset_shuffle = dataset[rand_inx,]
train.inx = 1:ceiling(.85 * nrow(dataset_shuffle))
    # val.inx = ceiling(.8 * nrow(dataset_shuffle)): ceiling(.9 * nrow(dataset_shuffle))
val.inx = ceiling(.85 * nrow(dataset_shuffle)):nrow(dataset_shuffle)

write.table(file=paste0(output.dir, "/dataset_train.txt"),x = dataset_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/dataset_val.txt"),x = dataset_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
