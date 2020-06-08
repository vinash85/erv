## convert to Seurat object and preprocess
library(data.table)
library(magrittr)
library(ggplo2)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(cowplot)
library(harmony)
library(uwot)
library(parallel)
library(EnhancedVolcano)
library(tidyverse)






## version 2 


myRowVar <- function(x) {
 (proxyC::rowSds(x))^2
}


mypreprocess_heterogeneous = function(X) {

    gene_var = myRowVar(X) 
    gene_var = ifelse(gene_var > 0, gene_var, NA)
    df = data.frame(gene = rownames(X), gene_mean = rowMeans(X),
                    zero_proportion = proxyC::rowZeros(X)/ncol(X),
                    gene_var = gene_var)
    df$samplesize = ncol(X)
    df = compute_test_statistic(df)
    return(df)
}

  

myone_level_clustering = function(subX, z_threshold, subsample= NULL) {
  require(bigpca)
  subdf = mypreprocess_heterogeneous(subX)
  # browser()
  # subdf = compute_test_statistic(subdf)
  features = subdf[subdf$zvalue > z_threshold, ]
  nullfeatures = data.frame(matrix(ncol = 11, nrow = 0))
  colnames(nullfeatures) = c("gene", "gene_mean", "zero_proportion",
                             "gene_var", "samplesize", "expected_pi", "se",
                             "minus_logp","zvalue", "subsetK", "K")
  if (nrow(features) < 10) {
    return(list(features = nullfeatures, pcs = NA, km = NA))
  }
  if (nrow(features) < 10) {
    return(list(features = nullfeatures, pcs = NA, km = NA,
                unscaled_pcs = NA,subdf = NA))
  }
     # browser()
     features.all = features
     # if(nrow(features) >100 )
     #  features = features[order(features$zvalue,decreasing=T)[1:100],] 

    # inx = which(rownames(subX) %in% features$gene)
    # mat = log(subX[inx, ] + 1)
    mat =  log1p(t(subX[features$gene,]))
    # mat = Seurat::LogNormalize(mat, scale.factor=1) %>% 
    #     scale(., center=T, scale=F)
    matMean = colMeans(mat)
    submat = mat
    
    if(!is.null(subsample)){
        size.curr = max(250, nrow(mat) * subsample)
        if(size.curr< nrow(mat))
            submat = mat[sample(nrow(mat), size=size.curr),]  
    }
  submat %<>% 
        scale(., center=T, scale=F)
  # pcs = NA
  # try(expr = {
  eigenvec = big.PCA(submat, pcs.to.keep =min(9, nrow(features) -1, ncol(submat) - 1), center = FALSE,  return.loadings = T)$loadings
  # browser()
  pcs = mat %*% eigenvec
  pcs.adjust = matMean %*% eigenvec
  pcs = sweep(pcs, 2, pcs.adjust, FUN = "-")
  # pcs = 
  # }, silent = TRUE)
  if (is.na(pcs[1])) {
      browser()
    return(list(features = nullfeatures, pcs = NA, km = NA,
                unscaled_pcs = NA,
                subdf = NA))
  } else {
    km = kmeans(pcs, 2, nstart = 10, iter.max = 50)
  }
  return(list(features = features.all,
              pcs = pcs,
              km = km,
              subdf = subdf))
}





myhippo = function(X, K = 20,
                 z_threshold = 2,
                 outlier_proportion = 0.001,
                 verbose = TRUE, subsample_pca = NULL) {

  if (outlier_proportion > 1 | outlier_proportion < 0) {
    stop("Outlier_proportion must be a number between 0 and 1.
         Default is 5%")
  }

  param = list(z_threshold = z_threshold,
               outlier_proportion = outlier_proportion,
               maxK = K)
  outlier_number = nrow(X) * outlier_proportion
  labelmatrix = matrix(NA, ncol(X), K)
  labelmatrix[, 1] = 1
  eachlevel = list()
  subX = X
  subXind = seq(ncol(X))
  withinss = rep(0, K)
  oldk = 1
  features = list()
  featuredata = list()
  for (k in 2:K) {
    thisk = myone_level_clustering(subX, z_threshold, subsample = subsample_pca)
    if (is.na(thisk$features$gene[1])) {
      if(verbose){
        message("not enough important features left; terminate the procedure")
      }
      labelmatrix = labelmatrix[, seq((k - 1))]
      break
    }
    if (nrow(thisk$features) < outlier_number) {
      if(verbose){
        message("not enough important features; terminate the procedure")
      }
      labelmatrix = labelmatrix[, seq((k - 1))]
      break
    }

    if (verbose) {message(paste0("K = ", k, ".."))}
    labelmatrix[, k] = labelmatrix[, k - 1]
    labelmatrix[subXind[thisk$km$cluster == 2], k] = k
    oneind = thisk$km$cluster == 1
    twoind = thisk$km$cluster == 2

    withinss[c(oldk,k)] = thisk$km$withins/thisk$km$size

    ind = which(table(thisk$km$cluster) <= 5)
    if (length(ind) >= 1){
      valid_indices = seq(k-1)[-ind]
      oldk = which(withinss == max(withinss[valid_indices]))
    }else{
      oldk = which.max(withinss[seq(k-1)])
    }
    if (sum(labelmatrix[, k] == oldk) < 2) {
      if(verbose){
        message("too few cells in one cluster; terminating the procedure")
      }
      labelmatrix = labelmatrix[, seq(k)]
      break
    }
    # subX = X[thisk$features$gene, which(labelmatrix[, k] == oldk)]
    subX = X[, which(labelmatrix[, k] == oldk)]
    subXind = which(labelmatrix[, k] == oldk)
    thisk$features$subsetK = oldk
    thisk$features$K = k
    features[[k - 1]] = thisk$features
  }
  hippo = list(X = X,
                                features = features,labelmatrix = labelmatrix,
                                z_threshold = z_threshold, param = param,
                                outlier_proportion = outlier_proportion)
  return(hippo)
}

setwd("~/project/deeplearning/icb/deepImmune/data_processing/ssgsea/imputation/simple")
source("~/liulab_home/softwares/HIPPO/R/hippo.R")


sce = readRDS("~/liulab_home/data/single_cell/Lee_data/lee.sce.RDS")
sce.cd3 = sce[,sce$CD3.status=="CD3+"]
sce.cd3.mat = assay(sce.cd3)
# rm(sce); gc()
# rm(sce.cd3); gc()

# temp  =sce.cd3.mat[1:1000,1:1000]
hippo.out = myhippo(
  sce.cd3.mat, 
  # temp, 
    K = 20, 
    z_threshold = 2, 
    outlier_proportion = 0.0001,
    verbose=TRUE, subsample_pca = .1)

hippo.out$cells = colnames(hippo.out$X) 
hippo.out.cd3pos = hippo.out 
hippo.out.cd3pos$X = NULL 
save(file="~/liulab_home/data/single_cell/Lee_data/cd3_plus.hippo.clustering.RData", hippo.out.cd3pos)


sce = readRDS("~/liulab_home/data/single_cell/Lee_data/lee.sce.RDS")
sce.cd3 = sce[,sce$CD3.status=="CD3-"]
sce.cd3.mat = assay(sce.cd3)
rm(sce); gc()
rm(sce.cd3); gc()
hippo.out.cd3neg = myhippo(sce.cd3.mat, 
    K = 30, 
    z_threshold = 2, 
    outlier_proportion = 0.0001,
    verbose=TRUE, subsample_pca = .1)

hippo.out.cd3neg$cells = colnames(hippo.out.cd3neg$X) 
hippo.out.cd3neg.cd3pos = hippo.out.cd3neg 
hippo.out.cd3neg$X = NULL 
save(file="~/liulab_home/data/single_cell/Lee_data/cd3_neg.hippo.clustering.RData", hippo.out.cd3neg)