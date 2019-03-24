
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

dataset.prefix = "Mel_PD1_Hugo"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_ICB_oxphos.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
dataset_ssgsea = sprintf("/liulab/asahu/data/ssgsea/xiaoman/expression/annot/%s.annot",dataset.prefix)
output.dir = sprintf("~/project/deeplearning/icb/data/%s", dataset.prefix)
pca_obj.RData = "/homes6/asahu/project/deeplearning/icb/data/tcga.blca/neoantigen/pca_obj.RData"
pca =TRUE
tpm=TRUE

library(data.table)
dir.create(output.dir, showWarnings = FALSE)

dataset_ssgsea = fread(dataset_ssgsea)
dataset_phenotype = fread(dataset_phenotype)
# genentech.env = local({load("/liulab/asahu/data/ssgsea/xiaoman/genentech.phenotype.RData");environment()})
    # phenotype_order = fread(phenotype_order)
dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
setnames(dataset_ssgsea, 1, "gene_name")
colnames(dataset_ssgsea_mat) = dataset_ssgsea$gene_name

dataset_phenotype$FMOne.mutation.burden.per.MB = -dataset_phenotype$Mutation
# identical(dataset_ssgsea$V1, pathway_order$pathway)
if(!tpm){
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$pathway]  
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$order]
    }else{
    pcgs = fread("/liulab/asahu/data/ssgsea/xiaoman/./pcg.txt")
           # dataset_ssgsea_mat = dataset_ssgsea_mat[,toupper(colnames(dataset_ssgsea_mat)) %in% toupper(pcgs$Gene)] 
    load("/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")
    # dataset_ssgsea_mat = dataset_ssgsea_mat[ ,common.genes]
  
    dataset_ssgsea_mat = impute.closest.gene(common.genes,dataset_ssgsea_mat)

    stopifnot(any(!is.na(dataset_ssgsea_mat)))

}


patient.name = rownames(dataset_ssgsea_mat)
patient.name = gsub(patient.name, pattern="-", replacement=".")
rownames(dataset_ssgsea_mat) = patient.name

# phenotype data
setnames(dataset_phenotype, 1, "patient.name")
dataset_phenotype$patient.name = gsub(dataset_phenotype$patient.name, pattern="-", replacement=".")
only_in_phenotype = setdiff(dataset_phenotype$patient.name, patient.name)
only_in_ssgsea = setdiff( patient.name, dataset_phenotype$patient.name)
common.patients = intersect(patient.name, dataset_phenotype$patient.name)
dataset_ssgsea_sel = dataset_ssgsea_mat[match(common.patients, patient.name), ] 

phenotype_sel = dataset_phenotype[match(common.patients, dataset_phenotype$patient.name)]
# phenotype_sel[1:2,]

colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern=" ", replacement="_")
colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern="-", replacement="_")
dataset_ssgsea_sel.back = dataset_ssgsea_sel





if(pca){
    load(pca_obj.RData)
            # pca_obj = NULL

    temp_out = get_pca(dataset_ssgsea_sel, pca_obj = pca_obj, scale=F) 
    pca_obj = temp_out$pca_obj
    pca_obj$len_selected = 50
    save(file=paste0(output.dir, "/pca_obj.RData"), pca_obj)
    pca_out_sel = temp_out$pca_out[,seq(pca_obj$len_selected)]
    dataset_ssgsea_sel = pca_out_sel 
}

response = phenotype_sel$Response
phenotype.mat = data.table(
    Best_CR = ifelse(response =="1", 1, 0),
    Best_PR = ifelse(response =="0.5", 1, 0),
    Best_PD = ifelse(response =="10", 1, 0),
    Best_SD = ifelse(response =="0", 1, 0), 
    Response = ifelse(response =="1"|response =="0.5" , 1, 0), 
    "Response_CR/PR" = ifelse(response =="1"|response =="0.5" , 1, 0), 
    "Response_SD/PD" = ifelse(response =="0", 1, 0)
    )

extra.genes.inx = c("TGFB1", "TGFBR2", "KLRC1") 
extra.genes.ez = c("7040", "7048", "3821") 
extra.genes = dataset_ssgsea_sel.back[, extra.genes.inx]
colnames(extra.genes) = extra.genes.ez
pheno.feat = phenotype_sel[,setdiff(5:41,39),with=F]
feat.merged = cbind(dataset_ssgsea_sel, extra.genes, pheno.feat)
pheno.merged = cbind(phenotype_sel[,which(colnames(phenotype_sel)!= "Response"), with=F],phenotype.mat)


genentech.pheno =  fread("~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/dataset_phenotype.txt")
genentech.feat =  fread("~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/dataset_ssgsea.txt")

feat.cols.na = setdiff(colnames(genentech.feat), colnames(feat.merged))
setdiff(colnames(genentech.pheno),colnames(pheno.merged))

feat.matched = as.matrix(feat.merged)
feat.matched = feat.matched[,match(colnames(genentech.feat),colnames(feat.matched))]
colnames(feat.matched) = colnames(genentech.feat)
feat.matched = apply(feat.matched,2,as.numeric)
feat.matched[,"SNV.Neoantigens"] = feat.matched[,"FMOne.mutation.burden.per.MB"]

pheno.matched = as.matrix(pheno.merged)
pheno.matched = pheno.matched[,match(colnames(genentech.pheno),colnames(pheno.matched))]
colnames(pheno.matched) = colnames(genentech.pheno)
pheno.matched = data.table(pheno.matched)
cols = setdiff(colnames(pheno.matched), "cancertype")
pheno.matched[ , (cols) := lapply(.SD, as.mynumeric), .SDcols = cols ]


write.table(file=paste0(output.dir, "/dataset_ssgsea.txt"),x = feat.matched,
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/dataset_phenotype.txt"),x = pheno.matched,
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/ssgsea_val.txt"),x = feat.matched,
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/phenotype_val.txt"),x = pheno.matched,
    row.names = F, col.names =T,  sep="\t", quote=F )




# imputation 


    feat.matched = fread("/homes6/asahu/project/deeplearning/icb/data/Mel_CTLA4_VanAllen/dataset_ssgsea.txt")
    pheno.matched = fread("/homes6/asahu/project/deeplearning/icb/data/Mel_CTLA4_VanAllen/dataset_phenotype.txt")
    predicted= fread("/homes6/asahu/project/deeplearning/icb//data/Mel_CTLA4_VanAllen/Neoantigen/val_prediction.csv",  skip=1)
    predicted.all = predicted[,-seq(ncol(predicted)/2), with=F]


    cols = c('SNV.Neoantigens',
        'Indel.Neoantigens',
        'Silent.Mutation.Rate',
        'Nonsilent.Mutation.Rate',
        'Number.of.Segments',
        'Fraction.Altered',
        'Aneuploidy.Scor',
        'HR')

    setnames(predicted.all, 1:8, cols)
    feat.removed = feat.matched[,setdiff(colnames(feat.matched),cols),with=F]
    pheno.removed = pheno.matched[,setdiff(colnames(pheno.matched),cols),with=F]
    Neoantigen.burden.per.MB = qnorm.array(feat.removed$Neoantigen.burden.per.MB)
    predicted.all$SNV.Neoantigens=ifelse(is.na(Neoantigen.burden.per.MB),predicted.all$SNV.Neoantigens,Neoantigen.burden.per.MB)
    FMOne.mutation.burden.per.MB = qnorm.array(feat.removed$FMOne.mutation.burden.per.MB)
    predicted.all$Silent.Mutation.Rate=ifelse(is.na(FMOne.mutation.burden.per.MB),predicted.all$Silent.Mutation.Rate,FMOne.mutation.burden.per.MB)

    feat.removed$Neoantigen.burden.per.MB = predicted.all$SNV.Neoantigens
    feat.removed$FMOne.mutation.burden.per.MB = predicted.all$Silent.Mutation.Rate




    feat.imputed = cbind(feat.removed, predicted.all)
    feat.imputed = feat.imputed[,colnames(feat.matched),with=F]
    pheno.imputed = cbind(pheno.removed, predicted.all)
    pheno.imputed = pheno.imputed[,colnames(pheno.matched),with=F]
    dir.create("/homes6/asahu/project/deeplearning/icb/data/Mel_CTLA4_VanAllen/imputed")
    write.table(file=paste0( "/homes6/asahu/project/deeplearning/icb/data/Mel_CTLA4_VanAllen/imputed/ssgsea_val.txt"),x = feat.imputed,
        row.names = F, col.names =T,  sep="\t", quote=F )
    write.table(file=paste0("/homes6/asahu/project/deeplearning/icb/data/Mel_CTLA4_VanAllen/imputed/phenotype_val.txt"),x = pheno.imputed,
        row.names = F, col.names =T,  sep="\t", quote=F )
