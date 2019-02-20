# standardise
minmax <- function(x) (x - min(x))/(max(x) - min(x))
x_train <- apply(ais[,1:11], 2, minmax)

# PCA
pca <- prcomp(x_train)

#
dataset_ssgsea_sel.old = dataset_ssgsea_sel
dataset_pca <- prcomp(scale(t(dataset_ssgsea_sel)),center = TRUE)
#Checking output of pca. prcomp function returns standard deviation (sdev), rotation and loadings.
aa = summary(dataset_pca)

dataset_ssgsea_sel = dataset_pca$rotation[,1:200]