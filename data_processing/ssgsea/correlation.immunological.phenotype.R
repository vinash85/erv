# please run process.dataset.imputation.tcga.R
seq(nrow(dataset_ssgsea_sel.back)
library(parallel)
out = apply(msi.neoantigen, 2, function(tt){
xx = mclapply(seq(ncol(dataset_ssgsea_sel.back)), function(uu) 
	c(cor.test(tt, dataset_ssgsea_sel.back[,uu])[c("estimate", "p.value")]), mc.cores=32)

temp = do.call(rbind,xx)
# rownames(temp) = colnames(dataset_ssgsea_sel.back)
}

 )
cors.mat  = out
genes = colnames(dataset_ssgsea_sel.back)
 save(file="/liulab/asahu/data/ssgsea/xiaoman/ICB_datasets/correlated_genes_with_ICB_biomarkers_p.value.RData", cors.mat, genes)


# correlation of everything



cors.mat = apply(dataset[,66:169,with=F], 2, function(tt){
xx = mclapply(seq(ncol(dataset_ssgsea_sel.back)), function(uu) 
	c(cor.test(tt, dataset_ssgsea_sel.back[,uu])[c("estimate", "p.value")]), mc.cores=32)

temp = do.call(rbind,xx)
rownames(temp) = colnames(dataset_ssgsea_sel.back)
temp
}

 )
genes = colnames(dataset_ssgsea_sel.back)
save(file="/liulab/asahu/data/ssgsea/xiaoman/ICB_datasets/correlated_genes_with_ICB_biomarkers_p.value.RData", cors.mat, genes)


## testing combat 

 mat = matrix(rnorm(30000), 1000, 30)
 mat = mat - min(mat)
 my.groups =factor(c(rep(0,15), rep(1,15)))
 inp.dt= mat 
 num.samp.thr = 5
 require(sva)
 # require(DESeq2)
 
 require(edgeR)
 require(limma)
 # inp.dt = vst(inp.dt)
 aa =  calcNormFactors(inp.dt, method = "TMM")
 inp.norm = inp.dt/aa
lcpm <- cpm(, log=TRUE)


    filter = apply(inp.dt, 1, function(x) length(x[x!=0])>=num.samp.thr)

    batch = as.numeric(my.groups)

    filtered = inp.dt[filter,]
    genes = rownames(filtered)[grep("^ENS", rownames(filtered))]
    dat0 = as.matrix(filtered)
    mod1 = model.matrix(~1, data=my.groups)
    # newdata = ComBat(log(filtered + .01), batch=batch, mod=mod1, par.prior = TRUE, prior.plots = FALSE)
    newdata = ComBat(filtered, batch=batch, mod=mod1, par.prior = TRUE, prior.plots = FALSE)
    inp.dt[filter,] = newdata
    rownames(newdata) = genes
    colnames(newdata) = common.samples
    newdata



## temp 

mat = normalize.expression(dataset_ssgsea_mat, num.samp.thr=0)
