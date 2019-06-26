
source("source.R")
dataset.prefix = "Getz_scRNA"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/getz/cell_label.csv"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
dataset_ssgsea =  "/liulab/asahu/data/ssgsea/xiaoman/getz/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"
output.dir = sprintf("~/project/deeplearning/icb/data/%s", dataset.prefix)
pca_obj.RData = "/homes6/asahu/project/deeplearning/icb/data/tcga.blca/neoantigen/pca_obj.RData"

pca_sel_obj.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/pca_sel_obj.RData"
ref.expression.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/ref.expression.RData"
ref.cancertype.RData = "~/project/deeplearning/icb/data/tcga/neoantigen.v2/ref.cancertype.RData"
library(data.table)
dir.create(output.dir, showWarnings = FALSE)
dataset_ssgsea_temp = fread(dataset_ssgsea, skip=2)
headers = fread(dataset_ssgsea, nrow=2)

dataset_phenotype = fread(dataset_phenotype)
xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea_temp[,seq(2,ncol(dataset_ssgsea_temp)),with=F]))
setnames(dataset_ssgsea_temp, 1, "gene_name")
colnames(dataset_ssgsea_mat) = dataset_ssgsea_temp$gene_name
rownames(dataset_ssgsea_mat) = unlist(headers[1])
dataset_ssgsea_mat = dataset_ssgsea_mat[-16292,]
tpm = T
pca = T
if(!tpm){
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$pathway]  
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$order]
    }else{
        pcgs = fread("/liulab/asahu/data/ssgsea/xiaoman/./pcg.txt")
           # dataset_ssgsea_mat = dataset_ssgsea_mat[,toupper(colnames(dataset_ssgsea_mat)) %in% toupper(pcgs$Gene)] 
        load("/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")
        dataset_ssgsea_mat1 = impute.closest.gene(common.genes,dataset_ssgsea_mat)
        dataset_ssgsea_mat = dataset_ssgsea_mat1[ ,common.genes] 

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
    ref.expression.cancertype = ref.expression 
    dataset_ssgsea_matched = match.expression.distribution(dataset_ssgsea_sel, ref.expression.cancertype)
    dataset_ssgsea_sel.back = dataset_ssgsea_matched

    dataset_ssgsea_sel = dataset_ssgsea_sel



    phenotype_sel.mod = phenotype_sel
    # generate cibersort results. 
    sc.assign = c("Bactivate", "Bnaive", "CD4Tmactivate", "CD4Tmresting", "CD8T", "Monocyte", "NKactivate", "Tfh", "Treg")
    cibersort.assign = c("Plasma_cells", "B_cells_naive", "T_cells_CD4_memory_activated", "T_cells_CD4_memory_resting", "T_cells_CD8", "Monocytes", "NK_cells_activated", "T_cells_follicular_helper", "T_cells_regulatory_(Tregs)")

    cell.types = phenotype_sel.mod$assign.ident 
    cibersort.map = cbind(sc.assign, cibersort.assign)
    cibersort.index = phenotype_order[5:27]
    cibersort.mat = matrix(0, nrow=nrow(phenotype_sel.mod), ncol=length(cibersort.index))
    colnames(cibersort.mat) = cibersort.index


    for (ii in seq(nrow(cibersort.map))) {
    map.curr = cibersort.map[ii,]
    cibersort.mat[cell.types==map.curr[1],map.curr[2]] = 1

    }
    cibersort.mat[,"Absolute_score"] = 1

    oxphos.dt = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/gad_oxphos_score.csv")
    oxphos.score = oxphos.dt[match(common.patients,object)]$Cor

    phenotype_sel.mod = cbind(phenotype_sel.mod, cibersort.mat,  dataset_ssgsea_sel[, c("CD274", "PDCD1")])
    phenotype_sel.mod = cbind(phenotype_sel.mod, oxphos.score)



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
cancertype[,"cancertype_SKCM"] = 1


dataset.small = cbind(phenotype.ext.mat[,1,with=F], general.pcs, pca_top, "MLH1" = dataset_ssgsea_matched[,"MLH1"], phenotype.ext.mat[,2:35,with=F], extra.genes, neoantigen.pheno, msi.pheno, cancertype)


dataset = dataset.small


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

save(file=paste0(output.dir, "/phenotype_data.RData"), phenotype_sel.mod)


## 
# read the model output
icb.phenotype = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/getz_val_prediction.csv")
# icb.samples = fread("~/project/deeplearning/icb/data/Getz_scRNA/sample_names_miao.txt")
icb.phenotype = icb.phenotype[unlist(icb.phenotype$sample_name) +1]

# response_information
resp = fread("GSE120575_patient_ID_single_cells.txt", skip=19)
xx = paste0("V",1:35)
colnames(resp)[1:35] = xx

# resp.patient
resp.matched=resp[match(phenotype_sel.mod$patient.name, resp$V2)]
response = resp.matched$V6
response.bin = ifelse(response=="Responder", 1, 0)
library(pROC)
cell.types = unique((phenotype_sel.mod$assign.ident) )
pre_post = resp$V5
pretreatment.samples = grep(pre_post, pattern="^Pre")
posttreatment.samples = grep(pre_post, pattern="^Post")

calc.stat = function(response.curr, value){
        tryCatch(
            auc(response.curr, value, levels=c(0,1)),
            error = function(e) NA
        )
}
plot.heatmap = function(dat, filename, height = 10, width =7){
  hc = hclust(as.dist(1-cor(dat, method="spearman", use="pairwise.complete.obs")), method="complete")
  hr = hclust(as.dist(1-cor(t(dat), method="spearman", use="pairwise.complete.obs")), method="complete")

  require(heatmap3)
  pdf( filename, height = height, width =width)

  heatmap3(dat, Rowv=as.dendrogram(hr),  Colv=as.dendrogram(hc), scale="none", balanceColor=T, showRowDendro=T ,   showColDendro=T, cexRow = .5, cexCol = 1)

  dev.off()
}

indexes = colnames(icb.phenotype)[2:115]

out = lapply(cell.types, function(cc) {
    inx = intersect(which(phenotype_sel.mod$assign.ident==cc), pretreatment.samples)
    response.curr= response.bin[inx]
    sapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat(response.curr, value.all[inx]) 
    }
        )
    })

out.dt = do.call(cbind, out)
colnames(out.dt) = cell.types
rownames(out.dt) = indexes
out.final = out.dt[rowSums(is.na(out.dt)) < 9,]
plot.heatmap(out.final, filename="pretreatment_aucs.pdf", height = 10, width =7)

# posttreatment
out = lapply(cell.types, function(cc) {
    inx = intersect(which(phenotype_sel.mod$assign.ident==cc), posttreatment.samples)
    response.curr= response.bin[inx]
    sapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat(response.curr, value.all[inx]) 
    }
        )
    })

out.dt = do.call(cbind, out)
colnames(out.dt) = cell.types
rownames(out.dt) = indexes
out.final = out.dt[rowSums(is.na(out.dt)) < 9,]
plot.heatmap(out.final, filename="posttreatment_aucs.pdf", height = 10, width =7)

## pre vs post in all
pre_or_post = rep(NA, nrow(phenotype_sel.mod))
pre_or_post[posttreatment.samples] = 1
pre_or_post[pretreatment.samples] = 0
out = lapply(cell.types, function(cc) {
    inx = which(phenotype_sel.mod$assign.ident==cc)
    response.curr= pre_or_post[inx]
    sapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat(response.curr, value.all[inx]) 
    }
        )
    })

out.dt = do.call(cbind, out)
colnames(out.dt) = cell.types
rownames(out.dt) = indexes
out.final = out.dt[rowSums(is.na(out.dt)) < 9,]
plot.heatmap(out.final, filename="post_vs_pretreatment_aucs_all.pdf", height = 10, width =7)

# responders 
out = lapply(cell.types, function(cc) {
    inx = intersect(which(phenotype_sel.mod$assign.ident==cc), which(response=="Responder"))
    response.curr= pre_or_post[inx] 
    sapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat(response.curr, value.all[inx]) 
    }
        )
    })

out.dt = do.call(cbind, out)
colnames(out.dt) = cell.types
rownames(out.dt) = indexes
out.final = out.dt[rowSums(is.na(out.dt)) < 9,]
plot.heatmap(out.final, filename="post_vs_pretreatment_aucs_responders.pdf", height = 10, width =7)
responders.pre_post = out.dt

# non responders
out = lapply(cell.types, function(cc) {
    inx = intersect(which(phenotype_sel.mod$assign.ident==cc), which(response=="Non-responder"))
    response.curr= pre_or_post[inx] 
    sapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat(response.curr, value.all[inx]) 
    }
        )
    })

out.dt = do.call(cbind, out)
colnames(out.dt) = cell.types
rownames(out.dt) = indexes
out.final = out.dt[rowSums(is.na(out.dt)) < 9,]
plot.heatmap(out.final, filename="post_vs_pretreatment_aucs_nonresponders.pdf", height = 10, width =7)
non_responders.pre_post = out.dt

nonresponders.pre_post = out.dt

## subtract responders vs non-responders
diff.responders.pre_post=  responders.pre_post - nonresponders.pre_post
plot.heatmap(diff.responders.pre_post, filename="post_vs_pretreatment_aucs_responders_vs_nonresponders.pdf", height = 10, width =7)

