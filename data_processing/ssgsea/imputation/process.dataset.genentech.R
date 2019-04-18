
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

options(error=recover)
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Genetech_expression_TPM.txt"

pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_ICB_oxphos.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/genentech.tpm/neoantigen.v2/"
pca_obj.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2//pca_obj.RData"
pca_sel_obj.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/pca_sel_obj.RData"
ref.expression.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/ref.expression.RData"
ref.cancertype.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/ref.cancertype.RData"
pca = T
tpm = T
library(data.table)
dir.create(output.dir, showWarnings = FALSE)

dataset_ssgsea = fread(dataset_ssgsea)
pathway_order = fread(pathway_order)
dataset_phenotype = fread(dataset_phenotype)
xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
genentech.env = local({load("/liulab/asahu/data/ssgsea/xiaoman/genentech.phenotype.RData");environment()})
    # phenotype_order = fread(phenotype_order)
dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
setnames(dataset_ssgsea, 1, "gene_name")
colnames(dataset_ssgsea_mat) = dataset_ssgsea$gene_name


# identical(dataset_ssgsea$V1, pathway_order$pathway)
tpm = T
if(!tpm)
{
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$pathway]  
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$order]
}else
{
    pcgs = fread("/liulab/asahu/data/ssgsea/xiaoman/./pcg.txt")
           # dataset_ssgsea_mat = dataset_ssgsea_mat[,toupper(colnames(dataset_ssgsea_mat)) %in% toupper(pcgs$Gene)] 
    load("/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")
    dataset_ssgsea_mat = dataset_ssgsea_mat[ ,common.genes] 
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

dataset_ssgsea_norm = normalize.expression(dataset_ssgsea_sel)
load(ref.expression.RData)
load(ref.cancertype.RData)
ref.expression.cancertype = ref.expression[ref.cancertype=="BLCA",] 
dataset_ssgsea_matched = match.expression.distribution(dataset_ssgsea_sel, ref.expression.cancertype)
dataset_ssgsea_sel.back = dataset_ssgsea_matched

dataset_ssgsea_sel = dataset_ssgsea_matched
            ######
            # imputing  response
            #########
# genetech.patients = intersect(grep(phenotype_sel$patient.name, pattern="^SAM", value=T), rownames(dataset_ssgsea_sel))
phenotype_sel.mod = phenotype_sel
phenotype_sel.mod[, Response:=as.double(Response)]
phenotype_sel.mod[is.na(Response) & (vital_status == 1) & (survive < 3)]$Response = 0
phenotype_sel.mod[is.na(Response) & (survive > 7)]$Response = 1
dataset_ssgsea_sel = dataset_ssgsea_sel
phenotype_order[length(phenotype_order)] = "Response" # last is response
phenotype_mat =  phenotype_sel.mod
temp = setdiff(phenotype_order, colnames(phenotype_mat))
temp.mat = matrix(NA, ncol=length(temp), nrow=nrow(phenotype_mat))
colnames(temp.mat) =temp
phenotype_mat = cbind(phenotype_mat, temp.mat)
phenotype.ext.mat = phenotype_mat[,match(phenotype_order, colnames(phenotype_mat)),with=F ]


if(pca){
    load(pca_obj.RData)

    temp_out = get_pca(dataset_ssgsea_sel, pca_obj = pca_obj, scale=F) 
    pca_obj = temp_out$pca_obj
    pca_obj$len_selected = 50
    save(file=paste0(output.dir, "/pca_obj.RData"), pca_obj)
    general.pcs = temp_out$pca_out[,seq(pca_obj$len_selected)]
    # dataset_ssgsea_sel = pca_out_sel 
}

reorder = match(rownames(dataset_ssgsea_sel), rownames(genentech.env$phenotype.feature.mat))
sels = c(1:8, 11:12)
f1 = genentech.env$phenotype.feature.mat[reorder,sels]
extra.genes.inx = c("TGFB1", "TGFBR2", "KLRC1")
f2 = t(genentech.env$cc[match(extra.genes.inx, genentech.env$bb$Symbol) ,])[reorder,]
colnames(f2)  = extra.genes.inx
pheno1 = genentech.env$response.mat[reorder,]

pheno.feat = phenotype.ext.mat[,5:34,drop=F]
new.feature = cbind(dataset_ssgsea_sel, f1, f2, pheno.feat)
new.pheno = cbind(phenotype.ext.mat, f1, f2, pheno1)
# dataset_ssgsea_sel = new.feature
# phenotype.ext.mat = new.pheno


# top pca
load(pca_sel_obj.RData) ## contains pca_sel_obj and top.genes 
temp_out = get_sel_pca(dataset_ssgsea_sel, top.genes, pca_sel_obj, scale=F)
pca_sel_obj = temp_out$pca_obj
pca_sel_obj$len_selected = 10
pca_top = temp_out$pca_out[,seq(pca_sel_obj$len_selected)]
colnames(pca_top) = paste0(colnames(pca_top), ".sel")


# extra genes
extra.genes.inx = c("TGFB1", "TGFBR2", "KLRC1") 
extra.genes.ez = c("7040", "7048", "3821") 
extra.genes = dataset_ssgsea_matched[, extra.genes.inx]
colnames(extra.genes) = extra.genes.inx

datasets.tcga = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2/dataset_train.txt")
# neoantigen.pheno
neoantigen.pheno.inx = colnames(datasets.tcga)[100:127]
neoantigen.pheno = matrix(NA, nrow=nrow(phenotype.ext.mat), ncol = length(neoantigen.pheno.inx))
colnames(neoantigen.pheno) = neoantigen.pheno.inx
reorder = match(rownames(dataset_ssgsea_sel), rownames(genentech.env$genentech.pheno))
neoantigen.pheno[,"Silent.Mutation.Rate"] = genentech.env$genentech.pheno$FMOne.mutation.burden.per.MB[reorder]
neoantigen.pheno[,"SNV.Neoantigens"] = genentech.env$genentech.pheno$Neoantigen.burden.per.MB[reorder]

# msi.pheno
msi.pheno.inx = colnames(datasets.tcga)[128:137]
msi.pheno = matrix(NA, nrow=nrow(phenotype.ext.mat), ncol = length(msi.pheno.inx))
colnames(msi.pheno) = msi.pheno.inx

#cancertype.pheno
msi.pheno.inx = colnames(datasets.tcga)[128:137]
msi.pheno = matrix(NA, nrow=nrow(phenotype.ext.mat), ncol = length(msi.pheno.inx))
colnames(msi.pheno) = msi.pheno.inx


cancertype.inx = colnames(datasets.tcga)[138:169]
cancertype = matrix(0, nrow=nrow(phenotype.ext.mat), ncol = length(cancertype.inx))
colnames(cancertype) = cancertype.inx
cancertype[,"cancertype_BLCA"] = 1


dataset.small = cbind(phenotype.ext.mat[,1,with=F], general.pcs, pca_top, "MLH1" = dataset_ssgsea_matched[,"MLH1"], phenotype.ext.mat[,2:35,with=F], extra.genes, neoantigen.pheno, msi.pheno, cancertype)

## adding genetech specific columns 
imputed.response = phenotype.ext.mat$Response
all.phenos = genentech.env$genentech.pheno[reorder,]

genetech.specific.cols = cbind(imputed.response, all.phenos)

dataset = cbind(dataset.small, genetech.specific.cols)

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




if(impute.neoantigen){
    pheno_inx = c("Neoantigen.burden.per.MB", "FMOne.mutation.burden.per.MB")
    reorder = match(rownames(dataset_ssgsea_sel), rownames(genentech.env$genentech.pheno))
                pheno1 = genentech.env$genentech.pheno[reorder,pheno_inx] ## changing to negative changed back to positive
                new.pheno = cbind(phenotype.ext.mat, pheno1)
                phenotype.ext.mat = new.pheno
            }
            if(imputed) source("data_processing/ssgsea/impute.genentech.neoantigen_plus.R")


            write.table(file=paste0(output.dir, "/dataset_ssgsea.txt"),x = dataset_ssgsea_sel,
                row.names = F, col.names =T,  sep="\t", quote=F )
            write.table(file=paste0(output.dir, "/dataset_phenotype.txt"),x = phenotype.ext.mat,
                row.names = F, col.names =T,  sep="\t", quote=F )



            rand_inx = sample(nrow(dataset_ssgsea_sel))
            dataset_ssgsea_sel_shuffle = dataset_ssgsea_sel[rand_inx,]
            phenotype.ext.mat_shuffle = phenotype.ext.mat[rand_inx,]
            train.inx = 1:ceiling(.8 * nrow(dataset_ssgsea_sel))
            val.inx = ceiling(.8 * nrow(dataset_ssgsea_sel)): ceiling( nrow(dataset_ssgsea_sel))

            library(pROC)

            aa = dataset_ssgsea_sel_shuffle
            bb= phenotype.ext.mat_shuffle
        # auc(bb$Response, aa$Neoantigen.burden.per.MB)
            auc(bb$Response, aa$SNV.Neoantigens)

        # test.inx = ceiling(.9 * nrow(dataset_ssgsea_sel)):nrow(dataset_ssgsea_sel)

            write.table(file=paste0(output.dir, "/ssgsea_train.txt"),x = dataset_ssgsea_sel_shuffle[train.inx,],
                row.names = F, col.names =T,  sep="\t", quote=F )
            write.table(file=paste0(output.dir, "/phenotype_train.txt"),x = phenotype.ext.mat_shuffle[train.inx,],
                row.names = F, col.names =T,  sep="\t", quote=F )

            write.table(file=paste0(output.dir, "/ssgsea_val.txt"),x = dataset_ssgsea_sel_shuffle[val.inx,],
                row.names = F, col.names =T,  sep="\t", quote=F )
            write.table(file=paste0(output.dir, "/phenotype_val.txt"),x = phenotype.ext.mat_shuffle[val.inx,],
                row.names = F, col.names =T,  sep="\t", quote=F )


