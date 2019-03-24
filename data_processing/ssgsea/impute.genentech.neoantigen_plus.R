# read data imputed from pre-trained model 
predicted= fread("../data/genentech.tpm/Neoantigen/val_prediction_pretrained.csv",  skip=1)
predicted.all = predicted[,-seq(ncol(predicted)/2), with=F]
actual = fread("../data/genentech.tpm/Neoantigen/dataset_phenotype.txt")


cols = c('SNV.Neoantigens',
'Indel.Neoantigens',
'Silent.Mutation.Rate',
'Nonsilent.Mutation.Rate',
'Number.of.Segments',
'Fraction.Altered',
'Aneuploidy.Scor',
'HR')
setnames(predicted.all, 1:8, cols)
# cor.test(actual$Neoantigen.burden.per.MB, predicted$V9, method="spearman")
Neoantigen.burden.per.MB = qnorm.array(actual$Neoantigen.burden.per.MB)
predicted.all$SNV.Neoantigens=ifelse(is.na(Neoantigen.burden.per.MB),predicted.all$SNV.Neoantigens,Neoantigen.burden.per.MB)
FMOne.mutation.burden.per.MB = qnorm.array(actual$FMOne.mutation.burden.per.MB)
predicted.all$Nonsilent.Mutation.Rate=ifelse(is.na(FMOne.mutation.burden.per.MB),predicted.all$Nonsilent.Mutation.Rate,FMOne.mutation.burden.per.MB)

phenotype.ext.mat = cbind(phenotype.ext.mat, predicted.all)
dataset_ssgsea_sel = cbind(dataset_ssgsea_sel, predicted.all)

