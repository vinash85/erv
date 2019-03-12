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
get_pca = function(data, pca_obj=NULL, center = T, scale = F){
	require(gmodels)
	if(is.null(pca_obj)){
			pca_obj <- big.prcomp(data,center = center, scale=scale)
		pca_out = pca_obj$x
		pca_obj$x = NULL ## to save space
	}else{
		aa = t(apply(data, 1, function(tt) {
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


