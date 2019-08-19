calc.stat.new = function(response.curr, value){
    aa = tryCatch(
        as.numeric(pROC::auc(response.curr, value, levels=c(0,1))),
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
# devtools::install_github("privefl/bigstatsr")

# parallel in R copied from https://privefl.github.io/blog/a-guide-to-parallelism-in-r/

calc.aucs.all = function(exp, inx,  response.curr, label = "gene", ngenes=NULL, ncluster=20){
    require(bigstatsr)
    calc.stat.new = function(response.curr, value){
        aa = tryCatch(
            as.numeric(pROC::auc(response.curr, value, levels=c(0,1))),
            error = function(e) NA
        )
        bb = tryCatch(
            wilcox.test(value[response.curr==0], value[response.curr==1], levels=c(0,1))$p.value,
            error = function(e) NA
        )
        c(aa,bb)
    }
    response.curr = FBM(nrow=1,ncol = length(inx), init=response.curr[inx])
    exp.curr = FBM(nrow=length(inx),  ncol = ncol(exp), init=exp[inx, ])
    out =  FBM(nrow=ncol(exp.curr),  ncol = 2, init=0) 
    library(foreach)
    library(doMC)
    registerDoMC(cores = ncluster)
    # cl <- parallel::makeCluster(ncluster)
    # doParallel::registerDoParallel(cl)
   foreach(tt = seq(ncol(exp.curr)), .combine = 'c') %dopar% {
        
        out[tt, ] = c(calc.stat.new(c(response.curr[1,]), c(exp.curr[,tt])))
        # out[tt,1] = aa[1]
        # out[tt,2] = aa[2]
        # aa[1]
        NULL

            }
    # parallel::stopCluster(cl)
    
    # out = mclapply(seq(ncol(exp.curr)), function(tt)  calc.stat.new(response.curr, exp.curr[,tt]), 
                   # mc.cores=32, mc.allow.recursive =T)
    # aucs.dt = do.call(rbind, out)
    aucs.dt = data.table(out[])
    aucs.dt$marker = colnames(exp)
    aucs.dt$label = label 
    aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    aucs.dt
}
# load("/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei.markers.RData")
# cell.type = chenfei.markers
load("/liulab/asahu/data/ssgsea/xiaoman/getz/response.bin.RData")
annot = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/SKCM_ICM_meta.txt")

load("/liulab/asahu/data/ssgsea/xiaoman/getz/icb.expression.matched.RData")
load("/liulab/asahu/data/ssgsea/xiaoman/getz/all.tcga.genes.RData")
phenotype_sel.mod = annot
setnames(phenotype_sel.mod, 1, "sample.name")
phenotype_sel.mod = phenotype_sel.mod[match(rownames(icb.expression.matched),sample.name)]

cell.types = phenotype_sel.mod[match(rownames(icb.expression.matched), sample.name)]$ assign.ident.2
resp = fread("/liulab/asahu/data/ssgsea/xiaoman/getz/GSE120575_patient_ID_single_cells.txt", skip=19)
xx = paste0("V",1:35)
colnames(resp)[1:35] = xx
resp$V2 = gsub(resp$V2, pattern="-", replacement=".")
resp.matched=resp[match(phenotype_sel.mod$sample.name, resp$V2)]
response = resp.matched$V6
response.bin = ifelse(response=="Responder", 1, 0)
#cell.types = unique((phenotype_sel.mod$assign.ident2) )
pre_post = resp$V5
pretreatment.samples = grep(pre_post, pattern="^Pre")
posttreatment.samples = grep(pre_post, pattern="^Post")

common.genes.scrna = intersect(all.tcga.genes, colnames(icb.expression.matched))
scrna.expr = icb.expression.matched[,common.genes.scrna]
# scrna.expr = icb.expression.matched[,1:10]
library(pROC)
library(foreach)
# registerDoMC(cores = 32)
response.curr = response.bin 
# out.all = foreach(cell.type = cell.types) %do% {
cell.types.unqiue = unique(cell.types) 
out.all = list()
for(cell.type in cell.types.unqiue[seq(8, length(cell.types.unqiue))]){ 
    print(cell.type) 
    inx = intersect(which(phenotype_sel.mod$assign.ident.2==cell.type), pretreatment.samples)
    topgene.Pretreatment = calc.aucs.all(scrna.expr, inx,  response.curr, label=sprintf("%s", cell.type))
    topgene.Pretreatment$treat = "Pre"
    inx = intersect(which(phenotype_sel.mod$assign.ident.2==cell.type), posttreatment.samples)
    topgene.Posttreatment = calc.aucs.all(scrna.expr, inx,  response.curr, label=sprintf("%s", cell.type))
    topgene.Posttreatment$treat = "Post"
    inx =which(phenotype_sel.mod$assign.ident.2==cell.type)
    topgene.Alltreatment = calc.aucs.all(scrna.expr, inx,  response.curr, label=sprintf("%s", cell.type))
    topgene.Alltreatment$treat = "All"
    out.all[[cell.type]] = rbind(topgene.Pretreatment, topgene.Posttreatment, topgene.Alltreatment)
    
}
topaucs.genes.list = out.all

topaucs.genes = do.call(rbind, out.all)
topaucs1 = topaucs.genes[aucs>0.79]

topaucs2 = do.call(rbind,lapply(out.all, function(tt) 
    rbind(
        tt[treat=="Pre"][order(aucs,decreasing=T)][V2<1E-7][1:10],
        tt[treat=="Post"][order(aucs,decreasing=T)][V2<1E-7][1:5],
        tt[treat=="All"][order(aucs,decreasing=T)][V2<1E-7][1:5])
    
))
topaucs2 = topaucs2[!is.na(V1)]
topaucs3 = topaucs.genes[aucs>0.75][V2<1E-10]
topaucs4 = topaucs.genes[V2<1E-15]
topaucs.final = rbind(topaucs1, topaucs2, topaucs3, topaucs4)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei1.marker.topaucs.genes.list.RData", topaucs.genes.list)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei1.marker.topaucs.final.RData", topaucs.final)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei1.marker.topaucs.genes.RData", topaucs.genes)
# save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/temp1.RData", out.all)
# xx = local({load("/liulab/asahu/data/ssgsea/xiaoman/getz/temp1.RData", ); environment()})
# yy = local({load("/liulab/asahu/data/ssgsea/xiaoman/getz/temp.RData", ); environment()})
# out.all = c(xx$out.all, yy$out.all, tempxx)
# save(file="~/project/deeplearning/icb/data/Getz_scRNA/phenotype_sel.mod.RData", phenotype_sel.mod)
