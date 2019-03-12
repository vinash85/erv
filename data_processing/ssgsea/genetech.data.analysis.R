### check in dessert tumors, responders vs. non-responder there is diff in exp. of NK2GA

library(IMvigor210CoreBiologies)
data(cds)
cc = counts(cds)
bb = fData(cds)
aa = pData(cds)
gene.inx = which(bb$symbol == "KLRC1")
gene.inx = which(bb$symbol == "TGFB1")
gene.inx = which(bb$entrez_id == 3133)
grp1 = cc[gene.inx,which(aa$"Immune phenotype"=="desert" & aa$binaryResponse == "SD/PD" )]
grp2 = cc[gene.inx, which(aa$"Immune phenotype"=="desert" & aa$binaryResponse == "CR/PR")]
grp1 = cc[gene.inx,which( aa$binaryResponse == "SD/PD" )]
grp2 = cc[gene.inx, which( aa$binaryResponse == "CR/PR")]
wilcox.test(grp1, grp2 )


## process phenotype data.
new.cols = colnames(aa)
new.cols = gsub(" ", new.cols, replacement=".")
new.cols[c(1,2,21, 22)]=c("Best", "Response", "OS", "Event")

genentech.pheno = aa
colnames(genentech.pheno) = new.cols

## select columns 
remove.cols = c(19, 20)
genentech.sel = genentech.pheno[, -remove.cols]
number.nas = apply(genentech.sel, 2,function(tt) sum(is.na(tt)) )
number.unique.values = apply(genentech.sel, 2,function(tt) length(unique(tt)) )
cbind(number.nas,number.unique.values, 1:23 )


response.phenotypes.inx = c("OS", "Event", "Response", "Best")
response.phenotypes = genentech.sel[,response.phenotypes.inx]
genentech.sel = genentech.sel[,setdiff(colnames(genentech.sel), response.phenotypes.inx)]
continuous.phenotypes.inx = c("FMOne.mutation.burden.per.MB", "Neoantigen.burden.per.MB")
continuous.phenotypes = genentech.sel[, continuous.phenotypes.inx]
genentech.sel = genentech.sel[,setdiff(colnames(genentech.sel), continuous.phenotypes.inx)]
genentech.sel = genentech.sel[,setdiff(colnames(genentech.sel), c("Sample.collected.pre-platinum", "Race"))]
# genentech.mat = genentech.sel
create_feature = function(xx, name, drop_na=T) {
	require(mltools)
	xx = as.factor(xx)
	if("MS2b2.1" %in% levels(xx)){
		new.lvls = levels(xx)
		new.lvls[grepl("MS2b2", new.lvls)] = "MS2b2" 
		levels(xx)  = new.lvls
	}
	lvls = table(xx) 
	common.lvls = names(lvls[lvls<=40])
	if(length(common.lvls) > 0 ){
		common.lab = paste(common.lvls,collapse="_")
		new.lvls = levels(xx)
		new.lvls[new.lvls%in%common.lvls] = common.lab
		levels(xx)  = new.lvls
	}

	lvls = sort(table(xx)) ## NA is inferred as most frequent phenotype
	temp = data.table(id=seq(length(xx)),xx)
	setnames(temp, 2, name)
	temp <- one_hot(temp)
	temp = as.matrix(temp[,-1,with=F])
	if(drop_na){
	temp[is.na(temp)] = 0 ## here is remove major class and automatically missing values are inferred as major class
	temp = temp[,-which.max(colSums(temp, na.rm=T)), drop=F] 
	}
	temp
}

temp = mclapply( seq(ncol(genentech.sel)), function(tt) create_feature(genentech.sel[,tt], colnames(genentech.sel)[tt]))

feature.mat = do.call(cbind, temp)

continuous.phenotypes.nan = as.matrix(continuous.phenotypes)
continuous.phenotypes.nan[is.na(continuous.phenotypes.nan)] = 0 

temp = mclapply( 3:4, function(tt) create_feature(response.phenotypes[,tt], colnames(response.phenotypes)[tt], drop_na=F))
response.mat = do.call(cbind, temp)

phenotype.feature.mat = cbind(feature.mat, continuous.phenotypes.nan)

save.image("/liulab/asahu/data/ssgsea/xiaoman/genentech.phenotype.RData")