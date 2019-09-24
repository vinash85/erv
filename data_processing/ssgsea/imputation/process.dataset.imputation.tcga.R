
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

dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/TCGA_ALLTPM.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/tcga_biom_oxphos.txt"
output.dir = "~/project/deeplearning/icb/data/tcga/scrna.v2/"
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




### extra phenotypes 


T_Cell_extra = c("IL10", "IDO", "TGFB1", "TGFB2", "TGFBR1", "TGFBR1", "CD37", "TLR", "Arginase")
APC_2 = c("A2AR", "VISTA", "B7_h3", "PDL1", "PDL2", "CD80", "CD86", "Galectin_9", "Ox40L", "CD40", "B7RP1", "CD70", "HVEM", "GITRL", "TNFSF9", "CD155", "CD112")
T_Cell_1 = c("CTLA4", "TIM3", "OX40", "CD40L", "ICOS", "CD27", "BTLA", "LAG3", "TCR", "KIR", "GITR", "TNFRSF9", "CD226", "TIGIT")
checkpoint.genes = unique(c(T_Cell_extra, APC_2, T_Cell_1))

all.tcga.genes = all.genes =dataset_ssgsea$gene_name
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/all.tcga.genes.RData", all.tcga.genes)
setdiff(checkpoint.genes, all.genes)
checkpoint.genes.1 = intersect(checkpoint.genes, all.genes)

#  [1] "IDO"        "TLR"        "Arginase"   "A2AR"       "VISTA"      "B7_h3"      "PDL1"       "PDL2"       "Galectin_9" "Ox40L"      "B7RP1"      "HVEM"       "GITRL"      "CD155"      "CD112"
# [16] "TIM3"       "OX40"       "CD40L"      "TCR"        "KIR"        "GITR"
checkpoint.genes.rescue = c("IDO1", "IDO2", "ARG1", "ARG2", "ADORA2A", "ADORA1", "VSIR", "CD276", "VTCN1", "JAK2", "STAT3", "CD80", "ICOSLG", "ICOS", "PVR", "CD226", "HAVCR2", "CD4", "PRF1", "FOXP3", "CD28", "LCK", "B2M")
pd1.genes = grep("^PDCD1",  all.genes, value=T)
Galectin_9.genes = grep("^LGALS",  all.genes, value=T)
jak.genes =c("JAK1", "JAK2", "JAK3")
stat.genes = grep("^STAT[0-9]",  all.genes, value=T)
TNF.genes  =c(grep("^TNFR",  all.genes, value=T), grep("^TNFS",  all.genes, value=T))
il2.gene = c("IL2", "PTPN2", grep("^IL2R",  all.genes, value=T))
il7.gene = grep("^IL7",  all.genes, value=T)
il4.gene = grep("^IL7",  all.genes, value=T)
il6.gene = grep("^IL6",  all.genes, value=T)
il10.gene = grep("^IL10",  all.genes, value=T)
HAVC.gene =grep("HAVC",  all.genes, value=T)
gzm.genes  =grep("^GZM",  all.genes, value=T)
traf.genes = grep("^TRAF",  all.genes, value=T)
nfk.genes = grep("^NFK",  all.genes, value=T)
cd40.genes = grep("^CD40",  all.genes, value=T)
igh.genes = grep("^IGH",  all.genes, value=T)
cd3.genes = grep("^CD3[A-Z]*$",  all.genes, value=T)
tra.genes = grep("^TR[A-B][C,D,V]",  all.genes, value=T)
kir.genes = grep("^KIR",  all.genes, value=T)
tgf.genes =grep("^TGF",  all.genes, value=T)
antigen.presentation.genes = grep("^HLA",  all.genes, value=T)
traf.genes = grep("^TR[A-B]F",  all.genes, value=T)
serpin.genes = grep("^SERPINB[1-9]$",  all.genes, value=T)
vegf.genes = grep("^VEGF",  all.genes, value=T)
tap.genes = c("TAP1", "TAP2", "TAPBP")

checkpoint.genes.semifinal = unique(c(checkpoint.genes.1, checkpoint.genes.rescue, pd1.genes, Galectin_9.genes, stat.genes, TNF.genes, il2.gene, il7.gene, il4.gene, il6.gene, il10.gene, HAVC.gene, gzm.genes, traf.genes, nfk.genes, cd40.genes, igh.genes, cd3.genes,tra.genes, kir.genes, tgf.genes, antigen.presentation.genes, traf.genes,serpin.genes, vegf.genes))
setdiff(checkpoint.genes.final, all.genes)
checkpoint.genes.final = intersect(checkpoint.genes.final, all.genes)



## scRNA data 

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

calc.aucs = function(exp, inx,  response.curr, label = "gene", ngenes=NULL){
 
    gene.select = genes.sel 
    response.curr = response.curr[inx]
    out = mclapply(gene.select, function(tt){
        value.all = exp[,tt]
        calc.stat.new(response.curr, value.all[inx]) 
        }, mc.cores=32, mc.allow.recursive =T
        )
    aucs.dt = do.call(rbind, out)
    aucs.dt = data.table(aucs.dt)
    aucs.dt$marker = dataset_ssgsea_temp$gene_name[gene.select]
    aucs.dt$label = label 
    aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    aucs.dt
}

calc.aucs.all = function(exp, inx,  response.curr, label = "gene", ngenes=NULL){
    response.curr = response.curr[inx]
    exp.curr = exp[inx, ]
    out = mclapply(seq(ncol(exp.curr)), function(tt)  calc.stat.new(response.curr, exp.curr[,tt]), 
        mc.cores=32, mc.allow.recursive =T
        )
    aucs.dt = do.call(rbind, out)
    aucs.dt = data.table(aucs.dt)
    aucs.dt$marker = colnames(exp)
    aucs.dt$label = label 
    aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    aucs.dt
}

common.genes.scrna = intersect(all.tcga.genes, colnames(icb.expression.matched))
scrna.expr = icb.expression.matched[,common.genes.scrna]
load("/liulab/asahu/data/ssgsea/xiaoman/getz/all.tcga.genes.RData")
require(doMC)
require(foreach)
registerDoMC(cores = 32)
response.curr = response.bin 
# out.all = foreach(cell.type = cell.types) %do% { 
out.all = list()
for(cell.type in unique(cell.types)){ 
    print(cell.type) 
    inx = intersect(which(phenotype_sel.mod$assign.ident==cell.type), pretreatment.samples)
    topgene.Pretreatment = calc.aucs.all(scrna.expr, inx,  response.curr, label=sprintf("%s", cell.type))
    topgene.Pretreatment$treat = "Pre"
    inx = intersect(which(phenotype_sel.mod$assign.ident==cell.type), posttreatment.samples)
    topgene.Posttreatment = calc.aucs.all(scrna.expr, inx,  response.curr, label=sprintf("%s", cell.type))
    topgene.Posttreatment$treat = "Post"
    inx =which(phenotype_sel.mod$assign.ident==cell.type)
    topgene.Alltreatment = calc.aucs.all(scrna.expr, inx,  response.curr, label=sprintf("%s", cell.type))
    topgene.Alltreatment$treat = "All"
    out.all[[cell.type]] = rbind(topgene.Pretreatment, topgene.Posttreatment, topgene.Alltreatment)

}
topaucs.genes.list = out.all

topaucs.genes = do.call(rbind, out.all)
topaucs1 = topaucs.genes[aucs>0.79]

topaucs2 = do.call(rbind,lapply(out.all, function(tt) 
    rbind(
        tt[treat=="Pre"][order(aucs,decreasing=T)][V2<1E-8][1:10],
        tt[treat=="Post"][order(aucs,decreasing=T)][V2<1E-8][1:10],
        tt[treat=="All"][order(aucs,decreasing=T)][V2<1E-8][1:5])

    ))
topaucs2 = topaucs2[!is.na(V1)]
topaucs3 = topaucs.genes[aucs>0.75][V2<1E-20]
topaucs4 = topaucs.genes[V2<1E-40]
topaucs.final = rbind(topaucs1, topaucs2, topaucs3, topaucs4)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/topaucs.genes.list.RData", topaucs.genes.list)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/topaucs.final.RData", topaucs.final)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/topaucs.genes.RData", topaucs.genes)



## all expression 
load("/liulab/asahu/data/ssgsea/xiaoman/getz/topaucs.final.RData")
immune.genes.compiled = intersect(all.genes, unique(c(topaucs.final$marker, checkpoint.genes.final)))

all.exp_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
patient.name = rownames(all.exp_mat)
patient.name = gsub(patient.name, pattern="-", replacement=".")
rownames(all.exp_mat) = patient.name
colnames(all.exp_mat) = all.genes
immune.genes.exps = all.exp_mat[common.patients, immune.genes.compiled]



## survival 
# "TCGA.V4.A9EW.01A" is in following format so need to convert 
common.patients.subs = substring(common.patients, 1,12)
library(readxl)
tcga_clinical = read_excel("/liulab/asahu/data/ssgsea/xiaoman/TCGA_gdc/TCGA-CDR-SupplementalTableS1.xlsx",sheet=1)
clinical.dt = data.table(tcga_clinical)
clinical.dt$patient.name = gsub(clinical.dt$bcr_patient_barcode, pattern="-", replacement=".")
clinical.dt = clinical.dt[match(common.patients.subs, patient.name)]

table3 = fread("/liulab/asahu/data/ssgsea/xiaoman/TCGA_gdc/TCGA-CDR_tabl3.csv")
setnames(table3, c(1,3,6,9,12), c('cancertype', paste0("sel.",c("OS", "PFI", "DFI", "DSS"))))
selstr="✓"

OS.types = table3[grep(table3$sel.OS, pattern=selstr)]$cancertype
clinical.dt[,OS.filtered:=ifelse(type %in% OS.types, OS, NA)]
PFI.types = table3[grep(table3$sel.PFI, pattern=selstr)]$cancertype
clinical.dt[,PFI.filtered:=ifelse(type %in% PFI.types, PFI, NA)]
DFI.types = table3[grep(table3$sel.DFI, pattern=selstr)]$cancertype
clinical.dt[,DFI.filtered:=ifelse(type %in% DFI.types, DFI, NA)]
DSS.types = table3[grep(table3$sel.DSS, pattern=selstr)]$cancertype
clinical.dt[,DSS.filtered:=ifelse(type %in% DSS.types, DSS, NA)]

clinical.final.dt = clinical.dt[, .(
    OS.time, OS.filtered, 
     PFI.time, PFI.filtered,
     DFI.time, DFI.filtered,
     DSS.time, DSS.filtered
    )]

# ERV genes 
# JCI121476.sdt12.txt 

ERV.tab = fread("/liulab/asahu/data/ssgsea/xiaoman/ERV_smith_jci/JCI121476.sdt12.txt")
common.patients.subs = substring(common.patients, 1,15)
ERV.patient.name = gsub(ERV.tab$Sample_ID, pattern="_", replacement=".")
ERV.mat=as.matrix(ERV.tab[,-1,with=F])
length(intersect(common.patients.subs, ERV.patient.name))
erv.pca.out = get_pca(ERV.mat, pca_obj = NULL, scale=F)

erv.pca_obj = erv.pca.out$pca_obj
erv.pca_obj$len_selected = 11
save(file=paste0(output.dir, "/erv.pca_obj.RData"), erv.pca_obj)
erv.pca_out_sel = erv.pca.out$pca_out[,seq(erv.pca_obj$len_selected)]
erv.pcs = erv.pca_out_sel
common.patients.erv = intersect(common.patients.subs, ERV.patient.name)
erv.pcs.common = erv.pcs[match(common.patients.erv,ERV.patient.name),]

erv.pcs.matched = matrix(NA, nrow=length(common.patients.subs), ncol = erv.pca_obj$len_selected)
erv.pcs.matched[match(common.patients.erv,common.patients.subs),] = erv.pcs.common


datasets.v2 = cbind(dataset, clinical.final.dt, erv.pcs.matched, immune.genes.exps)


write.table(file=paste0(output.dir, "/samples_name.txt"),x = common.patients,
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/dataset.txt"),x = datasets.v2,
    row.names = F, col.names =T,  sep="\t", quote=F )
# write.table(file=paste0(output.dir, "/dataset_phenotype.txt"),x = phenotype.ext.mat,
    # row.names = F, col.names =T,  sep="\t", quote=F )
rand_inx = sample(nrow(datasets.v2))
dataset_shuffle = datasets.v2[rand_inx,]
train.inx = 1:ceiling(.85 * nrow(dataset_shuffle))
    # val.inx = ceiling(.8 * nrow(dataset_shuffle)): ceiling(.9 * nrow(dataset_shuffle))
val.inx = ceiling(.85 * nrow(dataset_shuffle)):nrow(dataset_shuffle)

write.table(file=paste0(output.dir, "/dataset_train.txt"),x = dataset_shuffle[train.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )
write.table(file=paste0(output.dir, "/dataset_val.txt"),x = dataset_shuffle[val.inx,],
    row.names = F, col.names =T,  sep="\t", quote=F )

