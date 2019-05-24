
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
source("~/project/deeplearning/icb/deepImmune/source.R")

dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/TCGA_ssgsva.txt"
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/TCGA_ALLTPM.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/tcga_biom_oxphos.txt"
output.dir = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
# pca_obj.RData = "/homes6/asahu/project/deeplearning/icb/data/tcga.blca/neoantigen/pca_obj.RData"
tpm =T
pca = T

library(data.table)
dir.create(output.dir, showWarnings = FALSE)

dataset_ssgsea = fread(dataset_ssgsea)
pathway_order = fread(pathway_order)
xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
dataset_phenotype = fread(dataset_phenotype)
expression_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
tpm = T
if(!tpm){
    colnames(expression_mat) = dataset_ssgsea$V1
    expression_mat = expression_mat[,pathway_order$pathway]  
    expression_mat = expression_mat[,pathway_order$order]
    }else{
        colnames(expression_mat) = dataset_ssgsea$gene_name
        pcgs = fread("/liulab/asahu/data/ssgsea/xiaoman/./pcg.txt")
           # expression_mat = expression_mat[,toupper(colnames(expression_mat)) %in% toupper(pcgs$Gene)] 
        load("/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")
        expression_mat = expression_mat[ ,common.genes] 
        stopifnot(any(!is.na(expression_mat)))
    }


##BRCA

if(F){
    dataset_phenotype = dataset_phenotype[cancertype=="BLCA"]
}
patient.name = rownames(expression_mat)
patient.name = gsub(patient.name, pattern="-", replacement=".")
# patient.name = substring(patient.name, 1, 12)
rownames(expression_mat) = patient.name
## there are duplicates in patient names because same patient have multiple expression. 

# phenotype data
setnames(dataset_phenotype, 2, "patient.name")
dataset_phenotype$patient.name = gsub(dataset_phenotype$patient.name, pattern="-", replacement=".")
as.mynumeric = function(xx) as.numeric(ifelse(xx == '[Not Available]', NA, xx))
# dataset_phenotype1 =dataset_phenotype
cols = setdiff(colnames(dataset_phenotype)[-(1:2)], "cancertype")
dataset_phenotype[ , (cols) := lapply(.SD, as.mynumeric), .SDcols = cols ]
only_in_phenotype = setdiff(dataset_phenotype$patient.name, patient.name)
only_in_ssgsea = setdiff( patient.name, dataset_phenotype$patient.name)
common.patients = intersect(patient.name, dataset_phenotype$patient.name)
dataset_ssgsea_sel = expression_mat[match(common.patients, patient.name), ] 

phenotype_sel = dataset_phenotype[match(common.patients, dataset_phenotype$patient.name)]
# phenotype_sel[1:2,]

colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern=" ", replacement="_")
colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern="-", replacement="_")

normalize.survival = FALSE
prefix = unique(phenotype_sel$cancertype) 
if(normalize.survival){
    for (pre in seq(length(prefix))) {
        pre.curr = prefix[pre]
        inx.curr = which(phenotype_sel$cancertype == pre.curr)
        phenotype_sel$survive[inx.curr] = normalize.std(phenotype_sel$survive[inx.curr])


    }
}
phenotype_mat =  phenotype_sel[,-(1:2),with=F]
temp = setdiff(phenotype_order, colnames(phenotype_mat))
temp.mat = matrix(NA, ncol=length(temp), nrow=nrow(phenotype_mat))
colnames(temp.mat) =temp
phenotype.ext.mat = cbind(phenotype_mat, temp.mat)
phenotype.ext.mat = phenotype.ext.mat[,match(phenotype_order, colnames(phenotype.ext.mat)),with=F ]
dataset_ssgsea_sel = normalize.expression(dataset_ssgsea_sel, num.samp.thr =0) 
ref.expression = dataset_ssgsea_sel
ref.cancertype = phenotype.ext.mat$cancertype
save(file=paste0(output.dir, "/ref.expression.RData"), ref.expression)
save(file=paste0(output.dir, "/ref.cancertype.RData"), ref.cancertype)
dataset_ssgsea_sel.back = dataset_ssgsea_sel




if(pca){
    # load(pca_obj.RData)
    temp_out = get_pca(dataset_ssgsea_sel, pca_obj = NULL, scale=F, subsample=.4) 
        # temp_out = get_pca(dataset_ssgsea_sel, subsample=.2) 
    pca_obj = temp_out$pca_obj
    pca_obj$len_selected = 50
    save(file=paste0(output.dir, "/pca_obj.RData"), pca_obj)
    pca_out_sel = temp_out$pca_out[,seq(pca_obj$len_selected)]
    general.pcs = pca_out_sel
    # dataset_ssgsea_sel = pca_out_sel 
}
## add cancer type to the phenotye 


    barcode = substring(rownames(dataset_ssgsea_sel), 1,12)
        # library(readxl)
        # phenotype = excel_read("liulab/asahu/data/ssgsea/xiaoman/mmc2.xlsx")
    ## neoantigen 
    phenotype = fread("/liulab/asahu/data/ssgsea/xiaoman/mmc2.txt")
    cols = colnames(phenotype)
    cols[duplicated(cols)] = paste0(cols[duplicated(cols)], ".1")
    cols = gsub(cols, pattern=" ", replacement=".")
    colnames(phenotype) = cols

    phenotype$TCGA.Participant.Barcode = gsub(phenotype$TCGA.Participant.Barcode, pattern="-", replacement=".")
    sum(barcode %in% phenotype$TCGA.Participant.Barcode )
    reorder = match(barcode, phenotype$TCGA.Participant.Barcode)
    # feature.sel = c(3:36) ## if overall survival and progression free survival is considered 
    feature.sel = c(5:32)
    neoantigen.pheno = phenotype[reorder, feature.sel, with=F]

    ## msi 
    library(readxl)
    phenotype = read_excel("/liulab/asahu/data/ssgsea/xiaoman/msi.xlsx")
    phenotype = data.table(phenotype)
    phenotype$Barcode = gsub(phenotype$Barcode, pattern="-", replacement=".")
    sum(phenotype$Barcode %in% barcode)
    reorder = match(barcode, phenotype$Barcode)
    cols = colnames(phenotype)
    cols[duplicated(cols)] = paste0(cols[duplicated(cols)], ".1")
    cols = gsub(cols, pattern=" ", replacement=".")
    colnames(phenotype) = cols
    feature.sel = c(3:12)
    msi.pheno = phenotype[reorder, feature.sel, with=F]

    length(sort(table(phenotype.ext.mat$cancertype)))
    cancertype.dt = data.table(ID = rownames(dataset_ssgsea_sel), cancertype = as.factor(phenotype.ext.mat$cancertype))
    cancertype.pheno = mltools::one_hot(cancertype.dt, sparsifyNAs =T)

    cancertype.pheno = cancertype.pheno[,-1, with=F]

    extra.genes.inx = c("TGFB1", "TGFBR2", "KLRC1") 
    extra.genes.ez = c("7040", "7048", "3821") 
    extra.genes = dataset_ssgsea_sel.back[, extra.genes.inx]
    colnames(extra.genes) = extra.genes.inx

  

## select 10 genes per phenotype and take PC
msi.neoantigen = cbind(neoantigen.pheno, msi.pheno[,c("Total_nb_MSI_events", "MSI_exonic"), with=F]) 

cors = cor(msi.neoantigen, dataset_ssgsea_sel.back,  use = "pairwise.complete.obs")
top.cors = lapply(seq(nrow(cors)), function(tt) {
    xx = cors[tt,] 
    xx[order( abs(xx), decreasing = T)[1:20]] 
    }
)
names(top.cors)  = rownames(cors)
save(file="/liulab/asahu/data/ssgsea/xiaoman/ICB_datasets/top_correlated_genes_with_ICB_biomarkers.RData", top.cors)

 save(file="/liulab/asahu/data/ssgsea/xiaoman/ICB_datasets/correlated_genes_with_ICB_biomarkers.RData", cors)

top.genes = unique(unlist(lapply(top.cors, names)))

top.genes.extra = c("PMS2", "MSH6", "EPCAM", "MSH2")
top.genes = unique(c(top.genes, top.genes.extra))
setdiff(top.genes, common.genes)
top.expression = dataset_ssgsea_sel.back[,top.genes]

temp_out = get_sel_pca(top.expression, top.genes, scale=F)

pca_sel_obj = temp_out$pca_obj
pca_sel_obj$len_selected = 10
save(file=paste0(output.dir, "/pca_sel_obj.RData"), pca_sel_obj, top.genes)
pca_top = temp_out$pca_out[,seq(pca_sel_obj$len_selected)]

# cor.pcs = cor(msi.neoantigen, temp_out$pca_oat,  use = "pairwise.complete.obs")
# intersect(names(top.cors[["Indel.Neoantigens"]]), names(top.cors[["Total_nb_MSI_events"]]))
# pca_obj = temp_out$pca_obj
# pca_obj$len_selected = 10
colnames(pca_top) = paste0(colnames(pca_top), ".sel")
dataset = cbind(phenotype.ext.mat[,1,with=F], general.pcs, pca_top, "MLH1" = dataset_ssgsea_sel.back[,"MLH1"], phenotype.ext.mat[,2:35,with=F], extra.genes, neoantigen.pheno, msi.pheno, cancertype.pheno)




write.table(file=paste0(output.dir, "/samples_name.txt"),x = common.patients,
    row.names = F, col.names =T,  sep="\t", quote=F )
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



if(FALSE){
    phenotype_new_order = c("cancertype", phenotype_order[1:33], "oxphos_score",phenotype_order[34:41] )
    phenotype_order = phenotype_new_order
    save(file="/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData", phenotype_order)

}