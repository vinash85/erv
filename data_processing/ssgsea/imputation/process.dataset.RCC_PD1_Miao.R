
dataset.prefix = "RCC_PD1_Miao"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_ICB_oxphos.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
dataset_ssgsea = sprintf("/liulab/asahu/data/ssgsea/xiaoman/expression/annot/%s.annot",dataset.prefix)
output.dir = sprintf("~/project/deeplearning/icb/data/%s", dataset.prefix)
pca_obj.RData = "/homes6/asahu/project/deeplearning/icb/data/tcga.blca/neoantigen/pca_obj.RData"

pca_sel_obj.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/pca_sel_obj.RData"
ref.expression.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/ref.expression.RData"
ref.cancertype.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/ref.cancertype.RData"
library(data.table)
dir.create(output.dir, showWarnings = FALSE)

dataset_ssgsea = fread(dataset_ssgsea)
dataset_phenotype = fread(dataset_phenotype)
xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
setnames(dataset_ssgsea, 1, "gene_name")
colnames(dataset_ssgsea_mat) = dataset_ssgsea$gene_name


tpm = T
pca = T
if(!tpm){
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$pathway]  
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$order]
    }else{
        pcgs = fread("/liulab/asahu/data/ssgsea/xiaoman/./pcg.txt")
           # dataset_ssgsea_mat = dataset_ssgsea_mat[,toupper(colnames(dataset_ssgsea_mat)) %in% toupper(pcgs$Gene)] 
        load("/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")
        dataset_ssgsea_mat = impute.closest.gene(common.genes,dataset_ssgsea_mat)
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

    colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern=" ", replacement="_")
    colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern="-", replacement="_")
    dataset_ssgsea_norm = normalize.expression(dataset_ssgsea_sel)
    load(ref.expression.RData)
    load(ref.cancertype.RData)
    ref.expression.cancertype = ref.expression[ref.cancertype=="KIRC",] 
    dataset_ssgsea_matched = match.expression.distribution(dataset_ssgsea_sel, ref.expression.cancertype)
    dataset_ssgsea_sel.back = dataset_ssgsea_matched

    dataset_ssgsea_sel = dataset_ssgsea_sel

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
            # pca_obj = NULL

        temp_out = get_pca(dataset_ssgsea_sel, pca_obj = pca_obj, scale=F) 
        pca_obj = temp_out$pca_obj
        pca_obj$len_selected = 50
        save(file=paste0(output.dir, "/pca_obj.RData"), pca_obj)
        general.pcs = temp_out$pca_out[,seq(pca_obj$len_selected)]
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
    colnames(extra.genes) = extra.genes.inx


# top pca
load(pca_sel_obj.RData) ## contains pca_sel_obj and top.genes 
temp_out = get_sel_pca(dataset_ssgsea_sel.back, top.genes, pca_sel_obj, scale=F)
pca_sel_obj = temp_out$pca_obj
pca_sel_obj$len_selected = 10
pca_top = temp_out$pca_out[,seq(pca_sel_obj$len_selected)]
colnames(pca_top) = paste0(colnames(pca_top), ".sel")

datasets.tcga = fread("~/project/deeplearning/icb/data/tcga/neoantigen.v2/dataset_train.txt")
# neoantigen.pheno
neoantigen.pheno.inx = colnames(datasets.tcga)[100:127]
neoantigen.pheno = matrix(NA, nrow=nrow(phenotype.ext.mat), ncol = length(neoantigen.pheno.inx))
colnames(neoantigen.pheno) = neoantigen.pheno.inx

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
cancertype[,"cancertype_KIRC"] = 1


dataset.small = cbind(phenotype.ext.mat[,1,with=F], general.pcs, pca_top, "MLH1" = dataset_ssgsea_matched[,"MLH1"], phenotype.ext.mat[,2:35,with=F], extra.genes, neoantigen.pheno, msi.pheno, cancertype)

response.dt  = data.table(
    imputed.response = ifelse(response > 0,1,0),
    SD = ifelse(response ==0,1,0),
    PR = ifelse(response ==0.5,1,0),
    CR = ifelse(response ==1,1,0))
dataset = cbind(dataset.small, response.dt)


write.table(file=paste0(output.dir, "/dataset.txt"),x = dataset,
    row.names = F, col.names =T,  sep="\t", quote=F )

write.table(file=paste0(output.dir, "/sample_names.txt"),x = common.patients,
    row.names = F, col.names =T,  sep="\t", quote=F )
rand_inx = sample(nrow(dataset))
dataset_shuffle = dataset[rand_inx,]
train.inx = 1:ceiling(.85 * nrow(dataset_shuffle))
val.inx = ceiling(.85 * nrow(dataset_shuffle)):nrow(dataset_shuffle)

write.table(file=paste0(output.dir, "/dataset_train.txt"),x = dataset_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/dataset_val.txt"),x = dataset_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )



