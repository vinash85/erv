
source("~/project/deeplearning/icb/deepImmune/source.R")

dataset_ssgsea =  "/liulab/asahu/data/ssgsea/xiaoman/getz/dca/mean_norm.tsv"
dataset_phenotype =  "~/project/deeplearning/icb/data/Getz_scRNA/data/cell_label.csv"
options(error=recover)
dataset.prefix = "dca"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/getz/cell_label.csv"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))

load(phenotype_)

output.dir = sprintf("~/project/deeplearning/icb/data/Getz_scRNA/%s", dataset.prefix)
pca_obj.RData = "~/project/deeplearning/icb/data/tcga/scrna.v1//pca_obj.RData"

pca_sel_obj.RData = "~/project/deeplearning/icb/data/tcga/scrna.v1/pca_sel_obj.RData"
ref.expression.RData = "~/project/deeplearning/icb/data/tcga/scrna.v1/ref.expression.RData"
ref.cancertype.RData = "~/project/deeplearning/icb/data/tcga/scrna.v1/ref.cancertype.RData"
match.dataset.header = "~/project/deeplearning/icb/data/tcga/scrna.v1/dataset.txt"
library(data.table)
dir.create(output.dir, showWarnings = FALSE)


dataset_ssgsea_temp = fread(dataset_ssgsea)
dataset_phenotype = fread(dataset_phenotype)
# xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea_temp[,seq(2,ncol(dataset_ssgsea_temp)),with=F]))
setnames(dataset_ssgsea_temp, 1, "gene_name")
colnames(dataset_ssgsea_mat) = dataset_ssgsea_temp$gene_name
rownames(dataset_ssgsea_mat)[16292]
dataset_ssgsea_mat = dataset_ssgsea_mat[-16292,] # last column is NAs
rownames(dataset_ssgsea_mat) = colnames(dataset_ssgsea_temp)[-1]
# rownames(dataset_ssgsea_mat) = unlist(headers[1])[-1]
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
    phenotype_sel.mod = phenotype_sel
    dataset_ssgsea_norm = normalize.expression(dataset_ssgsea_sel)
    load(ref.expression.RData)
    load(ref.cancertype.RData)
    ref.expression.cancertype = ref.expression 
    dataset_ssgsea_matched = match.expression.distribution(dataset_ssgsea_sel, ref.expression.cancertype)
    dataset_ssgsea_sel.back = dataset_ssgsea_matched

    dataset_ssgsea_sel = dataset_ssgsea_sel.back



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
    phenotype_sel.mod = cbind(phenotype_sel.mod, oxphos_score=oxphos.score)



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

### add additonal columns to the dataset
match.dataset.header = colnames(fread(match.dataset.header, nrows=1))

# genes
checkpoint.genes = setdiff(intersect(match.dataset.header, dataset_ssgsea_temp$gene_name), colnames(dataset.small))
missing.genes = setdiff(match.dataset.header[seq(189,518)], dataset_ssgsea_temp$gene_name)


dataset_ssgsea_temp1 = dataset_ssgsea_temp[gene_name %in% checkpoint.genes]
gene_name_ord = dataset_ssgsea_temp$gene_name[dataset_ssgsea_temp$gene_name%in% checkpoint.genes]
dataset_ssgsea_mat1= t(as.matrix(dataset_ssgsea_temp1[,seq(2,ncol(dataset_ssgsea_temp1)),with=F]))
colnames(dataset_ssgsea_mat1) = gene_name_ord
rownames(dataset_ssgsea_mat1) = patient.name
dataset_checkpoint1 = dataset_ssgsea_mat1[common.patients,]
# phenotype data


dataset_genes_missing = impute.closest.gene(missing.genes,dataset_ssgsea_sel, exp2.dt=fread("/liulab/asahu/data/ssgsea/xiaoman/TCGA_ALLTPM.txt"))
dataset.curr = cbind(dataset, dataset_checkpoint1, dataset_genes_missing)
# NA columns
dataset.matched = match.matrix.col(dataset.curr, match.dataset.header)

dataset.new = data.table(dataset.matched)
write.table(file=paste0(output.dir, "/dataset.txt"),x = dataset.new,
    row.names = F, col.names =T,  sep="\t", quote=F )

write.table(file=paste0(output.dir, "/sample_names.txt"),x = common.patients,
    row.names = F, col.names =T,  sep="\t", quote=F )
rand_inx = sample(nrow(dataset.new))
dataset_shuffle = dataset.new[rand_inx,]
train.inx = 1:ceiling(.85 * nrow(dataset_shuffle))
val.inx = ceiling(.85 * nrow(dataset_shuffle)):nrow(dataset_shuffle)

write.table(file=paste0(output.dir, "/dataset_train.txt"),x = dataset_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/dataset_val.txt"),x = dataset_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )

save(file=paste0(output.dir, "/phenotype_data.RData"), phenotype_sel.mod)
save(file=paste0(output.dir, "/dataset_ssgsea_temp.RData"), dataset_ssgsea_temp) 
save(file=paste0(output.dir, "/headers.RData"), headers)
############### 
# read the model output
load("~/project/deeplearning/icb/data/Getz_scRNA/phenotype_data.RData")
icb.phenotype = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/getz_val_prediction.csv")
icb.phenotype = fread("~/project/deeplearning/icb/data/Getz_scRNA/val_prediction.csv")
# icb.samples = fread("~/project/deeplearning/icb/data/Getz_scRNA/sample_names_miao.txt")
icb.phenotype = icb.phenotype[unlist(icb.phenotype$sample_name) +1]
icb.phenotype.all = cbind(phenotype_sel.mod, icb.phenotype)
save(file="~/project/deeplearning/icb/data/Getz_scRNA/icb.phenotype.all.RData",icb.phenotype.all)
# response_information
resp = fread("~/project/deeplearning/icb/data/Getz_scRNA/data/GSE120575_patient_ID_single_cells.txt.gz", skip=19)
resp = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/GSE120575_patient_ID_single_cells.txt", skip=19)
xx = paste0("V",1:35)
colnames(resp)[1:35] = xx
resp$V2 = gsub(resp$V2, pattern="-", replacement=".")
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

  heatmap3(dat, Rowv=as.dendrogram(hr),  Colv=as.dendrogram(hc), scale="none", balanceColor=T, showRowDendro=F ,   showColDendro=F)
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
###############
## pre vs post in all
###############
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

###############
## subtract responders vs non-responders
###############
diff.responders.pre_post=  responders.pre_post - nonresponders.pre_post
diff.responders.pre_post = diff.responders.pre_post[rowSums(is.na(diff.responders.pre_post)) < 9,]
plot.heatmap(diff.responders.pre_post, filename="post_vs_pretreatment_aucs_responders_vs_nonresponders.pdf", height = 10, width =7)

###############
## evaluate expression of all genes in terms of auc and p-value of t-test 
###############

calc.stat.new = function(response.curr, value){
        aa = tryCatch(
            as.numeric(auc(response.curr, value, levels=c(0,1))),
            error = function(e) NA
        )
        bb = tryCatch(
            wilcox.test(value[response.curr==0], value[response.curr==1], levels=c(0,1))$p.value,
            error = function(e) NA
        )
        c(aa,bb)
}

icb.expression = t(dataset_ssgsea_temp[,2:16292, with=F])
colnames(icb.expression) =  dataset_ssgsea_temp$gene_name
rownames(icb.expression) =  gsub(rownames(dataset_ssgsea_mat), pattern="-", replacement=".")
phenotype_sel.mod$patient.name = gsub(phenotype_sel.mod$patient.name, pattern="-", replacement=".")
length(intersect(rownames(icb.expression), phenotype_sel.mod$patient.name))
icb.expression.matched = icb.expression[match(phenotype_sel.mod$patient.name, rownames(icb.expression)),]
temp1 = icb.expression.matched[,1:30]
colnames(icb.expression.matched) = dataset_ssgsea_temp$gene_name



plot.aucs.hist = function(inx, filename, title){
    cor.monocytes = cor(icb.expression.matched[inx,] ,response.bin[inx])
    genes.sel = order(abs(cor.monocytes),decreasing=T)[1:500]
    gene.select = unique( c(genes.sel, sample.int(length(cor.monocytes), 1000)))
# aa = auc( response.bin[inx], icb.expression.matched[inx, "HLA-G"]  )
    response.curr = response.bin[inx]
    out = mclapply(gene.select, function(tt){
        value.all = icb.expression.matched[,tt]
        calc.stat.new(response.curr, value.all[inx]) 
        }, mc.cores=32
        )

###############
# for each cell type create a figure for comparison of expression 
###############
    aucs.dt = do.call(rbind, out)
    aucs.dt = data.table(aucs.dt)
    aucs.dt$marker = dataset_ssgsea_temp$gene_name[gene.select]
    aucs.dt$label = "gene" 
    aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    out = mclapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat.new(response.curr, value.all[inx]) 
        }, mc.cores=32
        )
    di.aucs.dt = do.call(rbind, out)
    di.aucs.dt = data.table(di.aucs.dt)
    di.aucs.dt$marker = indexes
    di.aucs.dt$label = "signature"
    di.aucs.dt = di.aucs.dt[!is.na(V1)]
    di.aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    pre_treatment.aucs = rbind(aucs.dt, di.aucs.dt)
    pre_treatment.aucs = pre_treatment.aucs[order(aucs)]
    setnames(pre_treatment.aucs, "V2", "P")
    pre_treatment.aucs[,logP:=-log10(P)]
    require(ggrepel)
    m1 = di.aucs.dt[which(aucs > 0.7)]
    if(nrow(m1) > 25) m1 = di.aucs.dt[order(aucs,decreasing=T)][1:25]
    if(nrow(m1) < 2) m1 = di.aucs.dt[order(aucs,decreasing=T)[1:5]]
    m2 = aucs.dt[which(aucs > 0.7)]
    if(nrow(m2) > 20) m2 = aucs.dt[order(aucs,decreasing=T)[1:20]]
    if(nrow(m2) < 2) m2 = aucs.dt[order(aucs,decreasing=T)[1:5]]


    pre_treatment_subset = pre_treatment.aucs[marker %in% c(m1$marker, m2$marker)]
    p = ggplot(pre_treatment.aucs, aes(x = aucs, y = logP)) +
    geom_point(aes(color=as.factor(label)), alpha=0.5) +

    theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
    labs(x="AUC", y="Significance", title=title)+
    geom_text_repel(
        data = pre_treatment_subset,
        aes(x = aucs, y = logP, label = marker),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
        )  

    ggsave(p, file=filename)
}

###############
# pre-treatment macrophage
###############

require(doMC)
require(foreach)
registerDoMC(cores = 32)
out = foreach(cell.type = cell.types) %dopar% { 
    print(cell.type) 
    inx = intersect(which(phenotype_sel.mod$assign.ident==cell.type), pretreatment.samples)
    plot.aucs.hist(inx, filename = sprintf("/liulab/asahu/data/ssgsea/xiaoman/getz/aucs/pretreatment_%s.pdf", cell.type), title= sprintf("Pretreatment %s ", cell.type))
    inx = intersect(which(phenotype_sel.mod$assign.ident==cell.type), posttreatment.samples)
    plot.aucs.hist(inx, filename = sprintf("/liulab/asahu/data/ssgsea/xiaoman/getz/aucs/posttreatment_%s.pdf", cell.type), title= sprintf("Posttreatment %s ", cell.type))
    inx =which(phenotype_sel.mod$assign.ident==cell.type)
    plot.aucs.hist(inx, filename = sprintf("/liulab/asahu/data/ssgsea/xiaoman/getz/aucs/all_%s.pdf", cell.type), title= sprintf("All %s ", cell.type))

}




save(file="~/project/deeplearning/icb/data/Getz_scRNA/gadz.all.RData", icb.expression.matched, icb.phenotype, phenotype_sel.mod)

##
###############
###################
# create tsnes 
###################
###############
library(Rtsne)
library(ggplot2)
plotSNE <- function(data= data_tsne.merge,col = c(2:114),color.col = "condition",
                    title="t-SNE",size=0.25,do.discrete=T, filename=NULL, perplexity=30, theta=0.5, pca = FALSE, max_iter=5000, num_threads=32, tsne =NULL){
  set.seed(9)
  require(ggplot2)

  if(is.null(tsne)) tsne <- Rtsne(as.matrix(data[,col]), check_duplicates = FALSE, 
                pca = pca, perplexity=perplexity, theta=theta, dims=2, max_iter = max_iter, num_threads = num_threads)
  
  d_tsne_1 = as.data.frame(tsne$Y)
  d_tsne_1=cbind(d_tsne_1,data[,color.col,drop=F])
  
  ## plotting the results without clustering
  p=ggplot(d_tsne_1, aes(x=V1, y=V2)) +
    geom_point(size=size,aes_string(color=color.col), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("tSNE_1") + ylab("tSNE_2") +
    ggtitle(label = title) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
  if (do.discrete) {
    p<- p+ scale_colour_brewer(palette = "Set2")
  }
  ##theme(legend.position = "none")
  if(!is.null(filename)) ggsave(file=filename, p)
  list(tsne, p)
}

znorm = function(xx) (xx - mean(xx,na.rm=T))/sd(xx,na.rm=T)

plotSNE.array <- function(data= data_tsne.merge,col = c(2:114), color.cols = "condition",
    title="t-SNE",size=0.25,do.discrete=T, filename=NULL, perplexity=30, theta=0.5, pca = FALSE, max_iter=5000, normalize=TRUE, num_threads=32){
  set.seed(9)

  tsne <- Rtsne(as.matrix(data[,col]), check_duplicates = FALSE, 
    pca = pca, perplexity=perplexity, theta=theta, dims=2, max_iter = max_iter, num_threads = num_threads)
  dt1 = as.data.frame(tsne$Y)
  ps = list()
  for (color.col in color.cols) {
      if(normalize) data[[color.col]] = znorm(data[[color.col]])
      d_tsne_1=cbind(dt1,col=data[[color.col]], shape=data$shape)
      


      title.curr = sprintf("%s_%s", title, color.col)

  ## plotting the results without clustering
      p=ggplot(d_tsne_1, aes(x=V1, y=V2)) +
      geom_point(size=size,aes(color=col, shape=as.factor(shape)), alpha=0.8) +
      scale_color_gradient2(low = "blue", mid = "white",
                            high = "red", space = "Lab" ) + 
      guides(colour=guide_legend(override.aes=list(size=2))) +
      xlab("tSNE_1") + ylab("tSNE_2") +
      ggtitle(label = title.curr) +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
      if (do.discrete) {
        p<- p+ scale_colour_brewer(palette = "Set2")
    }
  ##theme(legend.position = "none")
    if(!is.null(filename)) {
      filename.curr = sprintf("%s_%s.pdf", filename, color.col)

      ggsave(file=filename.curr, p)
      ps[[color.col]]  = p
  }
}
list(d_tsne_1, ps)
}


color.clusters.features <- function(data, cluster,  color.cols = "condition",
    title="t-SNE",size=0.25,do.discrete=T, filename=NULL, normalize=TRUE){
    require(viridis)

  dt1 = as.data.frame(cluster)
  colnames(dt1) = c("V1", "V2")
  ps = list()
  for (color.col in color.cols) {
      if(normalize) data[[color.col]] = znorm(data[[color.col]])
      d_cluster_1=cbind(dt1,col=data[[color.col]], shape=data$shape)
      title.curr = sprintf("%s_%s", title, color.col)
  ## plotting the results without clustering
      p=ggplot(d_cluster_1, aes(x=V1, y=V2)) +
      geom_point(size=size,aes(color=col, shape=as.factor(shape)), alpha=0.7) +
      # scale_color_gradient2(low = "blue", mid = "white",
                            # high = "red", space = "Lab" ) +
                          
      # guides(colour=guide_legend(override.aes=list(size=2))) +
      xlab("Dim1") + ylab("Dim2") +
      ggtitle(label = title.curr) +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
      if (do.discrete) {
        p<- p+ scale_colour_brewer(palette = "Set2")
    }else{
         # p <- p+  scale_color_viridis() 
         p <- p+ scale_color_gradientn(colours = heat.colors(20, alpha=0.7, rev=T))

    }
  ##theme(legend.position = "none")
    if(!is.null(filename)) {
      filename.curr = sprintf("%s_%s.pdf", filename, color.col)

      ggsave(file=filename.curr, p)
      ps[[color.col]]  = p
  }
}
ps
}

# phenotype_sel.mod.back = phenotype_sel.mod
tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/tsnes/Monocyte"
dir.create(tsne.dir)
loc=grep(pattern = ".output",colnames(icb.phenotype))
loc2=grep(pattern = "embedding",colnames(icb.phenotype))
d1=icb.phenotype[,c(loc,loc2),with=F]
d2 = resp.matched[,.(V2, V5, V6,V7)]
setnames(d2, 1:4,c("sample.name", "patient.name", "response", "treatment"))
d3 = phenotype_sel.mod[,.(assign.ident, CD274, PDCD1, oxphos.score)]
data_tsne.merge = as.data.frame(cbind(d3, d2, d1))

perplexities = c(5, 10, 20, 25, 30, 45, 50,100)
thetas = c(0.0, 0.5, 0.75)
pcas= c(FALSE, TRUE)
require(doMC)
require(foreach)
registerDoMC(cores = 32)
foreach(perplexity =  perplexities) %dopar% {
foreach(theta =  thetas) %dopar% {
    for(pca in pcas){
        locx=which(data_tsne.merge$assign.ident=="Monocyte")
        title = sprintf("perplexity:%s pca:%s", perplexity, pca)
        all.p = plotSNE(data = data_tsne.merge[locx,], col=9:121, size = 2,do.discrete=F, title=title,
            color.col = "response", perplexity = perplexity, theta = theta, pca=pca, filename=sprintf("%s/%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", perplexity, pca, theta), max_iter=2000)

    }
}

}
## preplexity = 45, PCA=F, theta = 0 

locx=which(data_tsne.merge$assign.ident=="Monocyte")
all.p = plotSNE(data = data_tsne.merge[locx,], col=9:121, size = 2,do.discrete=F, title="Monocyte Response",
    color.col = "response", perplexity = 45, theta = 0, pca=FALSE, 
    filename=sprintf("%s/%s_%s_final.pdf", tsne.dir, "Monocyte", "response"), max_iter=5000)


# locx=which(data_tsne.merge$assign.ident=="Monocyte")
# all.p = plotSNE(data = data_tsne.merge[locx,], col=9:121, size = 1,do.discrete=F, title=title,
#         color.col = "response", filename=sprintf("%s/%s_%s.pdf", tsne.dir, "Monocyte", "response"))


# plot immue factors
tsne.dir = "~/project/deeplearning/icb/data/Getz_scRNA/data/Monocytes/tsnes/"
dir.create(tsne.dir)
immune.factors = colnames(data_tsne.merge)[9:121]
# data_tsne.merge$stroke = ifelse(data_tsne.merge$response == "Responder",  4, 0)
data_tsne.merge$shape = data_tsne.merge$response
final.tsne.p = plotSNE.array(data = data_tsne.merge[locx,], col=9:121, size = 4, do.discrete=F,title="Monocyte Response",
    color.cols = immune.factors, perplexity = 45, theta = 0, pca=FALSE, max_iter=5000,
        , filename=sprintf("%s/%s", tsne.dir, "Final"))


# plot expression 
## plot tsne for different monocyte markers 
#CD163 abd CD34 for M2 macrophage
# science 
# Mono1: CD14 VCAN S100A8 S100A9 FCN1 ITGAM
# Mono2: LAIR2 ASAH1 APOBEC3A TSPAN14 LIPA ITGAM
# Mono3 (also shares signature with Mono1): G0S2 CXCR1 CXCR2 NAMPT NEAT1 AL137655 CSF3R CD14 VCAN S100A8 S100A9 FCN1 ITGAM
# Mono4 (also shares signature with Mono1): PRF1 GNLY CTSW FGFBP2 IL2RB GZMA CD14 VCAN S100A8 S100A9 FCN1 ITGAM
genes.new = c("NOS2", "CIITA", "IL12", "ARG1", "YM1", "IL10", "MCR1", "FIZZ1", "Ccl2", "Ccl5", "IL10", "CD81", "H2EB")
cytokines.associated = c("IL18", "IL12", "IL1", "TNF", "IL10", "IL12", "IL10", "TGFB", "CCL2", "CCL5")
stat.molecules = c("STAT1", "STAT2", "STAT3", "STAT6", "STAT1")
irf.genes = c("IRF5", "IRF4", "IRF3", "P65", "P50")
addition.genes = c("PI3KCA", "AKT1", "MTOR") 
science.genes = unique(toupper(c("CD163", "CD34", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGAM", "LAIR2", "ASAH1", "APOBEC3A", "TSPAN14", "LIPA", "ITGAM", "G0S2", "CXCR1", "CXCR2", "NAMPT", "NEAT1", "AL137655", "CSF3R", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGAM", "", "PRF1", "GNLY", "CTSW", "FGFBP2", "IL2RB", "GZMA", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGAM")))
all.genes = unique(toupper(c(genes.new, cytokines.associated, stat.molecules, irf.genes, addition.genes, science.genes)))
sc.genes = colnames(icb.expression.matched)
setdiff(all.genes, sc.genes)
# "IL12"     "YM1"      "MCR1"     "FIZZ1"    "H2EB"     "IL1"      "TGFB"     "P65"      "P50"      "PI3KCA"   "AL137655"
alternative.name = c("IL12A", "IL12B", "CHIA", "MC1R", "RETNLB", "NFKB1", "RELA", "PIK3CA", "FAM157A")
setdiff(alternative.name, sc.genes)
HLA.genes = grep(sc.genes, pattern="^HLA", value=T)
TGF.genes = grep(sc.genes,  pattern="^TGFB", value=T)
csf.genes = grep(sc.genes, pattern="^CSF", value=T)
genes.final = intersect(sc.genes, unique(c(all.genes, HLA.genes, TGF.genes, csf.genes)))

icb.expression.interesting = icb.expression.matched[,genes.final]
data_tsne.genes = cbind(data_tsne.merge, icb.expression.interesting)

tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/tsnes/Monocyte/genes.unorm"
dir.create(tsne.dir)
genes.p = plotSNE.array(data = data_tsne.genes[locx,], col=9:121, size = 4, do.discrete=F,title="Monocyte Response",
    color.cols = genes.final, perplexity = 45, theta = 0, pca=FALSE, max_iter=5000,
        , filename=sprintf("%s/%s", tsne.dir, "Final"), normalize=F)



tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/tsnes/Monocyte/genes.tsne"
dir.create(tsne.dir)
registerDoMC(cores = 32)
foreach(perplexity =  perplexities) %dopar% {
foreach(theta =  thetas) %dopar% {
    for(pca in pcas){
        locx=which(data_tsne.genes$assign.ident=="Monocyte")
        title = sprintf("perplexity:%s pca:%s", perplexity, pca)
        all.p = plotSNE(data = data_tsne.genes[locx,], col=124:223, size = 2,do.discrete=F, title=title,
            color.col = "response", perplexity = perplexity, theta = theta, pca=pca, filename=sprintf("%s/%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", perplexity, pca, theta), max_iter=1000)

    }
}

}


data_tsne.all.genes = cbind(data_tsne.merge, icb.expression.matched)
perplexities = c(10, 20, 30, 45, 50,100)
thetas = c(0.0)
pcas= c(TRUE)
tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/tsnes/Monocyte/genes.all.tsne"
dir.create(tsne.dir)
foreach(perplexity =  perplexities) %dopar% {
foreach(theta =  thetas) %dopar% {
    for(pca in pcas){
        locx=which(data_tsne.all.genes$assign.ident=="Monocyte")
        title = sprintf("perplexity:%s pca:%s", perplexity, pca)
        all.p = plotSNE(data = data_tsne.all.genes[locx,], col=124:55860, size = 2,do.discrete=F, title=title,
            color.col = "response", perplexity = perplexity, theta = theta, pca=pca, filename=sprintf("%s/%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", perplexity, pca, theta), max_iter=1000)

    }
}

}
## UMAP
library(uwot) 

plotUMAP <- function(data, col = c(2:114), color.col = "condition",
                    title="UMAP",size=0.25, do.discrete=T, filename=NULL, n_neighbors =  15, learning_rate = 1, init = "spectral", min_dist = .01, pca = NULL,  n_threads=32){
  set.seed(9)

  umap.model <- umap(as.matrix(data[,col]), 
                pca = pca, n_neighbors =  n_neighbors, learning_rate = learning_rate, init = init, min_dist = min_dist, n_threads = n_threads, ret_model=T)

  
  d_umap_1 = as.data.frame(umap.model$embedding)
  d_umap_1=cbind(d_umap_1,data[,color.col,drop=F])
  
  ## plotting the results without clustering
  p=ggplot(d_umap_1, aes(x=V1, y=V2)) +
    geom_point(size=size,aes_string(color=color.col), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("umap_1") + ylab("umap_2") +
    ggtitle(label = title) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
  if (do.discrete) {
    p<- p+ scale_colour_brewer(palette = "Set2")
  }
  ##theme(legend.position = "none")
  if(!is.null(filename)) ggsave(file=filename, p)
  list(umap.model, p)
}

min_dists = c(0.1, 0.25, 0.5, 0.8, 0.99)
n_neighbors = c(2, 5, 10, 15, 20, 50, 100)
pcas = c(NULL, 20)
learning_rates = c(.1,.01,1) 
tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/umaps/Monocyte"
dir.create(tsne.dir)
foreach(min_dist =  min_dists) %dopar% {
foreach(n_neighbor =  n_neighbors) %dopar% {
foreach(pca =  pcas) %dopar% {
foreach(learning_rate =  learning_rates) %dopar% {
    all.p = plotUMAP(data = data_tsne.merge[locx,], col=9:121, size = 2,do.discrete=F, title="Monocyte Response",
        color.col = "response",
         n_neighbors =  n_neighbor, learning_rate = learning_rate, init = "spectral", min_dist = min_dist, pca = pca, 
        filename=sprintf("%s/%s_%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", min_dist, n_neighbor, ifelse(is.null(pca), 0, pca), learning_rate))
}}}}


foreach(min_dist =  min_dists) %dopar% {
foreach(n_neighbor =  n_neighbors) %dopar% {
foreach(learning_rate =  learning_rates) %dopar% {
    all.p = plotUMAP(data = data_tsne.merge[locx,], col=9:121, size = 2,do.discrete=F, title="Monocyte Response",
        color.col = "response",
         n_neighbors =  n_neighbor, learning_rate = learning_rate, init = "spectral", min_dist = min_dist, pca = NULL, 
        filename=sprintf("%s/%s_%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", min_dist, n_neighbor, "NULL", learning_rate))
}}}



## check if there are dendritic cells in monocytes
pca=TRUE
perplexity=30
theta=0
tsne.dir = "~/project/deeplearning/icb/data/Getz_scRNA/data/Monocytes/tsnes/all.genes"
title = "Monocyte all genes"
locx=which(data_tsne.merge$assign.ident=="Monocyte")
all.genes.tsne = plotSNE(data = data_tsne.all.genes[locx,], col=124:55858, size = 2,do.discrete=F, title=title,
            color.col = "response", perplexity = perplexity, theta = theta, pca=pca, filename=sprintf("%s/final_%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", perplexity, pca, theta), max_iter=1000, num_threads=10)


n_neighbor = 15; learning_rate = 1; min_dist = 0.01; pca = 50

all.umap.tsne= plotUMAP(data = data_tsne.all.genes[locx,], col=124:55858, do.discrete=F, title="Monocyte Response",
        color.col = "response",
         n_neighbors =  n_neighbor, learning_rate = learning_rate, init = "spectral", min_dist = min_dist, pca = NULL, 
        filename=sprintf("%s/%s_%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", min_dist, n_neighbor, "NULL", learning_rate))



# tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/tsnes/Monocyte/genes.all.tsne/staining"
tsne.dir = "~/project/deeplearning/icb/data/Getz_scRNA/data/Monocytes/tsnes/all.genes/staining"
dir.create(tsne.dir)
dc.genes = c("CD141", "CLEC9A", "CD1C", "CD11C", "AXL", "SIGLEC6", "PDC", "CD14", "CD16", "TLR2", "TLR4", "TLR7", "TLR9")
macrophage.genes =c("CD14", "CD40", "CD11B", "CD64", "EMR1", "MAC1", "MAC3", "CD68")
science.genes = unique(toupper(c("CD163", "CD34", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGAM", "LAIR2", "ASAH1", "APOBEC3A", "TSPAN14", "LIPA", "ITGAM", "G0S2", "CXCR1", "CXCR2", "NAMPT", "NEAT1", "AL137655", "CSF3R", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGAM", "", "PRF1", "GNLY", "CTSW", "FGFBP2", "IL2RB", "GZMA", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGAM")))
feature.sel  = unique(c(dc.genes, macrophage.genes))
setdiff(feature.sel, sc.genes) # "CD141" "CD11C" "CD16"  "CD11B" "CD64"  "MAC1"  "MAC3"
alt.genes = c("THBD", "ITGAX", "FCGR3A",  "ITGAM", "FCGR1A", "LGALS3")
setdiff(alt.genes, sc.genes) #
features.final = intersect(sc.genes, c(feature.sel, alt.genes))

data_tsne.all.genes$shape = data_tsne.merge$response
xx = color.clusters.features( data=data_tsne.all.genes[locx,], cluster=all.umap.tsne[[1]]$embedding,  color.cols = features.final,
    title="t-SNE",size=2, filename=sprintf("%s/%s", tsne.dir, "gene"), normalize=F, do.discrete=F)





### perform kmean clustering
# centers = 3; iter.max=10;  nstart=20
perform.kmeans = function(data, col = c(1:2), title = "kmean", centers = 3, iter.max=10,  nstart=20, size=0.25, filename = NULL){
    set.seed(20)
    colnames(data)[col] = c("V1", "V2")
    dataCluster <- kmeans(data[, col,drop=T], centers= centers, iter.max=iter.max,  nstart=nstart)
    data$clust = as.factor(dataCluster$cluster)

    p=ggplot(data, aes(x=V1, y=V2)) +
    geom_point(size=size,aes(color=clust), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("Dim1") + ylab("Dim2") +
    ggtitle(label = title) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
    # p<- p+ scale_colour_brewer(palette = "Set2")
  ##theme(legend.position = "none")
  if(!is.null(filename)) ggsave(file=filename, p)
  dataCluster
}


centers.all = 2:20
iter.max = 100
nstart = 20 
final.tsne.map = final.tsne.p[[1]]
kmean.dir = "~/project/deeplearning/icb/data/Getz_scRNA/data/Monocytes/tsnes/kmeans/"
dir.create(kmean.dir)
out = foreach(centers =  centers.all) %dopar% {
    all.p = perform.kmeans(
        data = final.tsne.map, col = c(1:2), centers = centers, iter.max=iter.max,  nstart=nstart, size=2, title="Monocyte kmean",
        filename=sprintf("%s/%s_%s.pdf", kmean.dir, "Monocyte", centers))
}


### perform clustering 
library(Seurat)
library(scater)

filename = "temp.pdf"
data.tsne= final.tsne.p[[1]]
data.exp = t(icb.expression.matched)
data.exp <- normalize(data.exp) # not -working
# Defining clusters and markers:
library(scran)
data.DI = data_tsne.all.genes[locx,9:121]
snn.gr <- buildSNNGraph(t(data.DI))
data.tsne$clust = factor(igraph::cluster_walktrap(snn.gr)$membership)
p=ggplot(data.tsne, aes(x=V1, y=V2)) +
geom_point(size=size,aes(color=clust), alpha=0.8) +
guides(colour=guide_legend(override.aes=list(size=2))) +
xlab("Dim1") + ylab("Dim2") +
ggtitle(label = title) +
theme_light(base_size=20) +
theme(axis.text.x=element_blank(),
  axis.text.y=element_blank()) 
if(!is.null(filename)) ggsave(file=filename, p)

library(scater)
monocyte.exp  = data.exp[,locx]
sce = SingleCellExperiment(
    assays = list(counts = monocyte.exp))


ave.counts <- rowMeans(counts(sce))
keep <- rowMeans(counts(sce)) >= 0.2
sum(keep)
sce1 = sce[keep, ]
sce1 = normalize(sce1)

# ave.counts <- rowMeans(counts(sce))
# keep <- ave.counts >= 1
# sum(keep)

# numcells <- nexprs(sce, byrow=TRUE)
# alt.keep <- numcells >= 10
# sum(alt.keep)


# sce <- SCESet(countData=monocyte.exp)
 
my.clusters = as.numeric(data.tsne$clust)
clust.col <- rainbow(max(my.clusters))
markers <- findMarkers(sce1, data.tsne$clust)

top.markers = unique(unlist(sapply(markers, function(aa) rownames(aa[aa$Top<= 15,]))))
top.exprs = sce1[top.markers,,drop=FALSE]
top.exprs = exprs(top.exprs)
heat.vals <- top.exprs - rowMeans(top.exprs)
# t(apply(top.exprs, 1,  function(tt) tt/max(tt)))
heat.vals.reorder =heat.vals[,order(my.clusters)]
my.clusters1 = my.clusters[order(my.clusters)]


# pdf("temp1.pdf")
# heatmap.2(heat.vals.reorder, col=bluered, symbreak=TRUE,  trace='none', cexRow=0.6,
#     ColSideColors=clust.col[my.clusters1], Rowv = FALSE, Colv=FALSE, dendrogram='none') 
#     legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters)), pch=16)
# dev.off()

pdf("temp2.pdf")
heatmap3(heat.vals.reorder,   cexRow=0.6,
    ColSideColors=clust.col[my.clusters1], Rowv = T, Colv=NA, showColDendro=F, showRowDendro=F) 
    legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters1)), pch=16)
dev.off()


dt.tsne = data.table(data.tsne)
dt.summary = dt.tsne[,sum(shape=="Responder")/sum(shape!="Responder"), by=clust]
dt.summary = dt.summary[order(clust)]

dt.summary$clust = factor(dt.summary$clust, levels=seq(7,1))
p = ggplot(data=dt.summary, aes(x=clust, y=V1)) + 
geom_col(aes(fill=V1)) + coord_flip() + theme_minimal() 
ggsave(file="temp3.pdf", p)


### plot average level deepImmune predictions. 
monocyte.immune.factors  = t(data_tsne.all.genes[locx,9:57]) 
monocyte.immune.norm = t(apply(monocyte.immune.factors, 1, qnorm.array))

monocyte.immune.reorder =monocyte.immune.norm[,order(my.clusters)]
my.clusters1 = my.clusters[order(my.clusters)]

pdf("~/project/deeplearning/icb/data/Getz_scRNA/data/Monocytes/immune_factor_heatmap.pdf")
heatmap3(monocyte.immune.reorder,   cexRow=0.6,
    ColSideColors=clust.col[my.clusters1], Rowv = T, Colv=NA, showColDendro=F, showRowDendro=F) 
    # legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters1)), pch=16)
dev.off()







## only use tsne data 

tsne.dir = "/liulab/asahu/data/ssgsea/xiaoman/getz/tsnes/Monocyte_emb"
dir.create(tsne.dir)

perplexities = c(5, 10, 20, 25, 30, 45, 50,100)
thetas = c(0.0, 0.5, 0.75)
pcas= c(FALSE, TRUE)
require(doMC)
require(foreach)
registerDoMC(cores = 32)
foreach(perplexity =  perplexities) %dopar% {
foreach(theta =  thetas) %dopar% {
    for(pca in pcas){
        locx=which(data_tsne.merge$assign.ident=="Monocyte")
        title = sprintf("perplexity:%s pca:%s", perplexity, pca)
        all.p = plotSNE(data = data_tsne.merge[locx,], col=58:121, size = 2,do.discrete=F, title=title,
            color.col = "response", perplexity = perplexity, theta = theta, pca=pca, filename=sprintf("%s/%s_%s_%s_%s_%s.pdf", tsne.dir, "Monocyte", "response", perplexity, pca, theta), max_iter=2000)

    }
}

}



aa = dataset_ssgsea_temp[,2:5,with=F]
colSums(aa)
bb =rowSums(dataset_ssgsea_mat)
quantile(bb,probs=seq(0,1,.05))

dataset_ssgsea_temp[,2:5,with=F]

uu = sapply(list.rsem, function(tt) {
    tt1 = basename(tt)
    gsub(tt1, pattern ="_L001.rsem.genes.results", replacement="")
    })
pat.name = gsub(resp$V2, pattern="_L001$", replacement="")
xx = setdiff(pat.name, uu)
yy = setdiff(uu,pat.name)
intersect(pat.name, uu)
 sum(grepl(list.rsem,pattern="_L001.rsem.genes.results"))


aa = dataset_ssgsea_mat[grep( rownames(dataset_ssgsea_mat), pattern="H5_P2_M41"),]
# in kraken 

###################
# create approx read counts
###################

approx.tpm2count = function(TPM_mat, effective_length, scale=.00107){
    round(apply(TPM_mat, 1, function(tt) tt*effective_length)*scale)
}


# library(GEOquery)
# gpl <- getGEO("GPL18573", GSEMatrix=T)
# gpl = Table(gpl)


# gpl <- getGEO("GSE120575")
# gpl = Table(gpl)

# library("biomaRt")  
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# map <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),   mart= ensembl)
# gene.id = gsub(temp$gene_id, perl=T, pattern="\\.[0-9]+$", replacement="")
# common.genes = intersect(gene.id, map$ensembl_gene_id)
# id.inx = match(common.genes, gene.id)
scale=  510.0275 * 1E-6
temp = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/A1_P1_M39_L001.rsem.genes.results")
genecode.id = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/gencode.v19.annotation.gtf", skip=5, header=F)
map  =data.table(t(sapply(genecode.id$V9,function(tt) strsplit(tt, split=";| ")[[1]][c(5,14)])))

map[,V3:=gsub(V1, pattern="\"", replacement="")]
map[,V4:=gsub(V2, pattern="\"", replacement="")]
map = map[!duplicated(V3)]
tpm.mat = 2^dataset_ssgsea_mat -1 
setkey(map, V3)
temp$gene_name = map[temp$gene_id]$V4
sum(colnames(tpm.mat) %in% temp$gene_name) ## sample x genes are column 
effective_length = temp[match(colnames(tpm.mat), gene_name)]$effective_length
effective_length[effective_length==0] =1
xx = rowSums(tpm.mat)
count.mat = approx.tpm2count(tpm.mat, effective_length, scale=scale)
output.dir1 = "~/project/deeplearning/icb/data/Getz_scRNA/impute/"
dir.create(output.dir1)
write.table(file=paste0(output.dir1, "/count.mat.txt"), x = count.mat,
    row.names = T, col.names =T,  sep="\t", quote=F)

## corelate with genentech response data.
#1. individual gene
#2. individual gene + monocyte
#3. combined genes 
#4. combined gene + monocyte 


 # CD8A, PD1, CTLA4, and LAG3 