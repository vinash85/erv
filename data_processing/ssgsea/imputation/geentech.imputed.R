

feat = fread("~/project/deeplearning/icb/data/genentech.tpm/neoantigen.v2/dataset.txt")
predicted = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/tensorboardLog/survival20190422-081755/genentech_validation.csv")


# samples.names = rownames(dataset_ssgsea_sel)
samples.names = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/genentech.imputed/samples.names.txt")$x
cols = colnames(predicted)
predicted.all = predicted[,grep(cols, pattern=".output$"), with=F]
cols = colnames(predicted.all)
setnames(predicted.all, cols, gsub(cols, pattern=".output", replacement=""))
cols = colnames(predicted.all)
cols_to_replace = cols[13:50]


SNV.Neoantigens = match.distribution.zeros(feat$SNV.Neoantigens,predicted.all$SNV.Neoantigens) 
SNV.Neoantigens=ifelse(is.na(SNV.Neoantigens),predicted.all$SNV.Neoantigens,SNV.Neoantigens)

Nonsilent.Mutation.Rate = match.distribution.zeros(feat$Silent.Mutation.Rate,predicted.all$Nonsilent.Mutation.Rate)
Nonsilent.Mutation.Rate=ifelse(is.na(Nonsilent.Mutation.Rate),predicted.all$Nonsilent.Mutation.Rate,Nonsilent.Mutation.Rate)



feat.new = cbind(feat[,-cols_to_replace, with=F], predicted.all[,cols_to_replace,with=F]) 
feat.new = feat.new[,colnames(feat),with=F]
feat.new$SNV.Neoantigens = SNV.Neoantigens
feat.new$Nonsilent.Mutation.Rate = Nonsilent.Mutation.Rate


genentech.env = local({load("/liulab/asahu/data/ssgsea/xiaoman/genentech.phenotype.RData");environment()})
reorder = match(samples.names, rownames(genentech.env$phenotype.feature.mat))
sels = c(1:8, 11:12)
f1 = genentech.env$phenotype.feature.mat[reorder,sels]
pheno1 = genentech.env$response.mat[reorder,]


dataset = cbind(feat.new, f1, pheno1)

ref.dir = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention"
output.dir = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/genentech.imputed"

write.dataset(output.dir,dataset, ref.dir)



# tcga data imputation and analysis 

feat = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/dataset.txt")
predicted = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/tensorboardLog/survival20190422-081755/tcga_validation.csv")

# samples.names = rownames(dataset_ssgsea_sel)
cols = colnames(predicted)
predicted.all = predicted[,grep(cols, pattern=".output$"), with=F]
cols = colnames(predicted.all)
setnames(predicted.all, cols, gsub(cols, pattern=".output", replacement=""))
cols = colnames(predicted.all)
cols_to_replace = cols[-1]

feat.new =feat
impute_col = function(old, pred){
    old = match.distribution.zeros(old,pred) 
    ifelse(is.na(old),pred,old)

    } 
for (colx in cols_to_replace){
    feat.new[[colx]] = impute_col(feat.new[[colx]], predicted.all[[colx]])
}



dataset = cbind(feat.new, f1, pheno1)

ref.dir = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention"
output.dir = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/attention/genentech.imputed"

write.dataset(output.dir,dataset, ref.dir)


