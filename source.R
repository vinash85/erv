## all source utlitity function for R processing 


normalize.std = function(tt){
    (tt - min(tt, na.rm=T))/(max(tt,na.rm=T) - min(tt, na.rm=T))
}

# standardise
minmax <- function(x) (x - min(x))/(max(x) - min(x))

big.prcomp = function(data, center=TRUE, scale=FALSE){
	data1 = scale(data, center=center, scale=scale)
	require(bigpca)
	out = list()
	pca = big.PCA(data1,  pcs.to.keep = ncol(data1))
	out$center = center
	if(center) out$center = colMeans(data)
	out$scale = scale
	if(scale) out$scale = apply(data,2,sd)
	out$rotation = pca$PCs
	out$x = data1 %*% out$rotation
	# out$x = pca$loadings
	out
}
get_pca = function(data, pca_obj=NULL, center = T, scale = F, subsample=1){
	require(gmodels)
	sds =  apply(data, 2,  sd)
	zeros.sds = which(sds == 0)
	sds = ifelse(sds==0, 1E-3, sds)
	pca_out = NULL
	# calculate pca 
	if(is.null(pca_obj)){

		if(subsample < 1){
			data.sub = data[sample(nrow(data),size = subsample*nrow(data)),]
			pca_obj <- big.prcomp(data.sub,center = center, scale=scale)

			}else{
				pca_obj <- big.prcomp(data,center = center, scale=scale)
				pca_out = pca_obj$x
			}
		pca_obj$x = NULL ## to save space
		pca_obj$sds = sds
	}
	if(is.null(pca_out)) {
		aa = t(apply(data, 1, function(tt) {
			## matching sds 
			tt = tt * pca_obj$sds/sds ## matching the distribution
			if(pca_obj$center[1]) 
				tt = tt-pca_obj$center
			if(pca_obj$scale[1]) 
				tt = tt/pca_obj$scale
			tt
		}))
		pca_out = aa %*% pca_obj$rotation 
	}
	list(pca_out=pca_out, pca_obj=pca_obj)
}


range01.norm = function(mat){
	m1 = min(mat,na.rm=T)
	m2 = max(mat,na.rm=T)
	if(m1==m2) m2 = m2 +1
	(mat - m1)/(m2-m1)
}

	
qnorm.array <- function(mat)
{
	mat.back = mat
	mat = mat.back[!is.na(mat.back)]
    mat = rank(mat,  rank, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat
    mat.back
}

impute.closest.gene = function(common.genes, dataset_ssgsea_mat){
	  genes.imputed = setdiff(common.genes, colnames(dataset_ssgsea_mat))
	  if(length(genes.imputed) > 0) {
	  gene1 = colnames(dataset_ssgsea_mat)
	  exp2.dt = fread("/liulab/asahu/data/ssgsea/xiaoman/Genetech_expression_TPM.txt")
	  exp2 = t(as.matrix(exp2.dt[,seq(2,ncol(exp2.dt)),with=F]))
	  setnames(exp2.dt, 1, "gene_name")
	  colnames(exp2 ) = exp2.dt$gene_name
	  impute = exp2[,genes.imputed]
	  only.genes = intersect(gene1, exp2.dt$gene_name)
	  dataset_ssgsea_mat = dataset_ssgsea_mat[,only.genes]
	  exp.present = exp2[,only.genes]
	  cors = cor(impute, exp.present)
	  genes.inx = apply(cors,1, 
	  	function(tt) ifelse(sum(!is.na(tt)), which.max(tt), NA)
	  	)

	  imputed = dataset_ssgsea_mat[,genes.inx]
	  imputed[is.na(imputed)] = 0
	  colnames(imputed) = genes.imputed
	  merged = cbind(dataset_ssgsea_mat, imputed) 
	  dataset_ssgsea_mat = merged
	}
	dataset_ssgsea_mat[,common.genes]
}


create.km = function( times1, times2,  labels=list(Del="not-rescued", Norm="rescued"), file="temp.pdf",ppt=F){
    require(data.table)
    require(survival)
    require(ggplot2)
    dt = data.table( rbind(times1, times2)  )
    dt$label = labels$Del
    dt$label[1:nrow(times1)] = labels$Norm
    setnames(dt, 1:2, c("time", "status"))
    sr.surv <- survfit(Surv(time,status==0)~label, data = dt)
    dt$marker <- c(rep(0,nrow(times1)), rep(1,nrow(times2)))
    dt$status1 = ifelse(dt$status==0,1,0)
    max_x = 3*max(times1[,1])/4 
    outf3 = logRank(times1, times2)
    aa <- ggsurv(sr.surv)
    text = paste0(
        drug, "\n",
        "P=", formatC(outf3[1], format="E", digits=2), "\n",
        "AUC=", formatC(outf3[8] - outf3[7],digits=2)
        )
    if(ppt){
      aa <-aa +
      annotate( "text", x = max_x, y = .6, label = text) 
  }else{
      aa <-aa +
      annotate( "text", x = max_x, y = .6, label = text) +
      theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())  
  }
  # pdf(file)
  # print(aa)
  # dev.off()
  aa 
}

