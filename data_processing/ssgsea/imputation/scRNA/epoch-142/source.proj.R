recompile_r_avinash = function()
    recompile_r(lib.name ="avinash", path="~/shortcuts/avinash", lib = "/homes6/asahu/R/x86_64-pc-linux-gnu-library/3.6")
recompile_r = function(path, lib.name, reload=T, ...){
    # tryCatch( eval(parse(text=paste( "detach(package:", lib.name, ", unload = TRUE)", sep=""))), error = function(e) print("not loaded... "))
    # lib.namer <- paste("package:",lib.name,sep="")
    # tryCatch(library.dynam.unload( lib.name, system.file(package = lib.name)), error = function(e) print("not loaded... "))
    remove.packages(lib.name,...)
    install.opts <- "--no-lock "
    install.packages(path, repos=NULL, clean = TRUE,INSTALL_opts = install.opts, ...) 
    if(reload) require(lib.name, character.only = TRUE)
}
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

plot.aucs.hist = function(inx, indexes, filename, title, icb.phenotype=icb.phenotype, aucs.dt = NULL){
    
    require(ggthemes)
    response.curr = response.bin[inx]
    if(is.null(aucs.dt)){
        cor.monocytes = cor(icb.expression.matched[inx,] ,response.bin[inx])
        genes.sel = order(abs(cor.monocytes),decreasing=T)[1:500]
        gene.select = unique( c(genes.sel, sample.int(length(cor.monocytes), 1000)))
        # aa = auc( response.bin[inx], icb.expression.matched[inx, "HLA-G"]  )
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
        
        aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    }
    aucs.dt = aucs.dt[,.(V1, V2, marker, label, aucs)]
    aucs.dt$label = "gene"
    aucs.dt$alpha = 0.35
    
    out = mclapply(indexes, function(tt){
        value.all = icb.phenotype[[tt]]
        calc.stat.new(response.curr, value.all[inx]) 
    }, mc.cores=32
    )
    di.aucs.dt = do.call(rbind, out)
    di.aucs.dt = data.table(di.aucs.dt)
    di.aucs.dt$marker = gsub(indexes, pattern = ".output$", replacement = "")
    di.aucs.dt$label = "signature"
    di.aucs.dt = di.aucs.dt[!is.na(V1)]
    di.aucs.dt[,aucs:=ifelse( V1 < 0.5, 1-V1, V1)]
    di.aucs.dt$alpha = 0.8
    pre_treatment.aucs = rbind(aucs.dt, di.aucs.dt)
    pre_treatment.aucs = pre_treatment.aucs[order(aucs)]
    setnames(pre_treatment.aucs, "V2", "P")
    pre_treatment.aucs[,logP:=-log10(P)]
    require(ggrepel)
    m1 = di.aucs.dt[which(aucs > 0.7)]
    if(nrow(m1) > 20) m1 = di.aucs.dt[order(aucs,decreasing=T)][1:25]
    if(nrow(m1) < 2) m1 = di.aucs.dt[order(aucs,decreasing=T)[1:5]]
    m2 = aucs.dt[which(aucs > 0.7)]
    if(nrow(m2) > 20) m2 = aucs.dt[order(aucs,decreasing=T)[1:20]]
    if(nrow(m2) < 2) m2 = aucs.dt[order(aucs,decreasing=T)[1:5]]
    
    
    pre_treatment_subset = pre_treatment.aucs[marker %in% c(m1$marker, m2$marker)]
    p = ggplot(pre_treatment.aucs, aes(x = aucs, y = logP)) +
        geom_point(aes(color=as.factor(label), alpha = alpha)) +
        
        theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
        ggthemes::scale_colour_tableau() + 
        labs(x="AUC", y="Significance", title=title)+
        geom_text_repel(
            data = pre_treatment_subset,
            aes(x = aucs, y = logP, label = marker),
            size = 3,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )  
    
    ggsave(p, file=filename, width =7, height = 7)
    p
}





plot.all.aucs = function(cwd, icb.phenotype.col.dt=icb.phenotype.col.dt, icb.phenotype=icb.phenotype){
    
    load("/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei1.marker.topaucs.genes.list.RData")
    # save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei1.marker.topaucs.final.RData", topaucs.final)
    load("/liulab/asahu/data/ssgsea/xiaoman/getz/chenfei1.marker.topaucs.genes.RData")
    aa = topaucs.genes.list$Monocyte[treat=="Pre"]
    
    auc.dir = sprintf("%s/aucs",cwd)
    dir.create(auc.dir)
    
    pheno.inx.start = icb.phenotype.col.dt[V1=="B_cells_naive.output"]$V2
    pheno.inx.end = tail(grep(icb.phenotype.col.dt$V1, pattern="output$"),1)
    indexes.output = colnames(icb.phenotype)[pheno.inx.start:pheno.inx.end]
    survival.inx.start = icb.phenotype.col.dt[V1=="OS.time.output"]$V2
    indexes1 = grep(colnames(icb.phenotype), pattern = "embedding", value = T)
    indexes.survival = colnames(icb.phenotype)[seq(survival.inx.start, survival.inx.start+3)]
    selected.outputs = indexes = c( indexes1, indexes.output, indexes.survival)
    
    require(doMC)
    require(foreach)
    registerDoMC(cores = 32)
    Pre.p = Post.p = All.p = list()
    library(parallel)
    # out = foreach(cell.type = cell.types) %dopar% { 
    for(cell.type in cell.types) { 
        print(cell.type)
        
        cell.type.str = gsub(cell.type, pattern="/", replacement=".or.")
        cell.type.str = gsub(cell.type.str, pattern="-", replacement=".or.")
        inx = intersect(which(phenotype_sel.mod$assign.ident.2==cell.type), pretreatment.samples)
        # aucs.dt = topaucs.genes[treat=="Pre" & label==cell.type]
        aucs.dt = (topaucs.genes.list[[cell.type]])[treat=="Pre"][order(V2)][1:2000]
        Pre.p[[cell.type]] = plot.aucs.hist(inx, indexes = indexes, filename =sprintf("%s/pretreatment_%s.pdf", auc.dir, cell.type.str), title= sprintf("Pretreatment %s", cell.type.str), aucs.dt = aucs.dt)
        
        inx = intersect(which(phenotype_sel.mod$assign.ident.2==cell.type), posttreatment.samples)
        aucs.dt = topaucs.genes.list[[cell.type]][treat=="Post"][order(V2)][1:2000]
        Post.p[[cell.type]] = plot.aucs.hist(inx, indexes = indexes, filename = sprintf("%s/posttreatment_%s.pdf", auc.dir, cell.type.str), title= sprintf("Posttreatment %s  ", cell.type.str), aucs.dt = aucs.dt)
        
        
        inx =which(phenotype_sel.mod$assign.ident.2==cell.type)
        aucs.dt = topaucs.genes.list[[cell.type]][treat=="All"][order(V2)][1:2000]
        All.p[[cell.type]] = plot.aucs.hist(inx, indexes = indexes, filename = sprintf("%s/Alltreatment_%s.pdf", auc.dir, cell.type.str), title= sprintf("All %s  ", cell.type.str), aucs.dt = aucs.dt)
        
    }
    
    
    
    biauc.dir = sprintf("%s/biaucs",cwd)
    dir.create(biauc.dir)
    icb.phenotype.output = icb.phenotype[,selected.outputs, with=F]
    setnames(icb.phenotype.output, colnames(icb.phenotype.output), gsub(colnames(icb.phenotype.output), pattern=".output$", replacement="") )
    mat = cbind(icb.expression.matched, as.matrix(icb.phenotype.output))
    Pre.all.aucs = plot.biauc(ps=Pre.p,  mat=mat, dir = biauc.dir, cell.types = phenotype_sel.mod$assign.ident.2, indexes =pretreatment.samples)
    Pre.all.aucs[["auc.ps"]] = Pre.p
    Post.all.aucs = plot.biauc(ps=Post.p,  mat=mat, dir = biauc.dir, cell.types = phenotype_sel.mod$assign.ident.2, indexes =posttreatment.samples)
    Post.all.aucs[["auc.ps"]] = Post.p
    All.all.aucs = plot.biauc(ps=All.p,  mat=mat, dir = biauc.dir, cell.types = phenotype_sel.mod$assign.ident.2, indexes = seq(nrow(mat)))
    All.all.aucs[["auc.ps"]] = All.p
    
    
    AUCs = list(Pre.all.aucs, Post.all.aucs,All.all.aucs)
    save(file=sprintf("%s/AUCs.RData",biauc.dir), AUCs)
    
}


# ## Characterization of embedding 
# 1. Correlation 
# 2. High /low wilcox test
# a. Immune factors
# b. genes 


# aa = data.table(temp.corr$t[2,], temp.corr$cor[2,], colnames(temp.corr$cor))[order(abs(V2))]
# cor.curr = cor(x=embedding.matrix, y=pheno.matrix, method = "spearman")
agg.cor.volcano = function(cors, Ps, markers, filename=NULL, title="emb"){
    
    df.val = data.table(val=cors, P=Ps, marker = markers)
    # setnames(df.cors, "V2", "P")
    df.val[,logP:=-log10(P)]
    require(ggrepel)
    m2 = df.val[which(P < 1E-20 & abs(val) > 0.3)]
    if(nrow(m2) > 20) m2 = df.val[order(P)[1:5]]
    if(nrow(m2) < 2) m2 = df.val[order(P)[1:2]]
    thr = ifelse(sum(df.val$P < 1E-3) > 5, 1E-3, 1E-2)
    df.val[,Significant:=ifelse(P < thr, "Significant", "Not-significant")]
    df.val$title = title
    # df_subset = df.val[marker %in%  m2$marker]
    df.val[,repel:=ifelse(marker %in%  m2$marker, T, F)]
    df.val
}

agg.tvalue.volcano = function(cors, tvalues, markers, filename=NULL, title="emb"){
    
    df.val = data.table(val=cors, tvalue=tvalues, marker = markers)
    require(ggrepel)
    m2 = df.val[which(abs(val) > 0.3)]
    if(nrow(m2) > 5) m2 = df.val[order(tvalue,decreasing = T)[1:5]]
    if(nrow(m2) < 2) m2 = df.val[order(tvalue,decreasing = T)[1:2]]
    thr = sort(df.val$tvalue, decreasing = T )[10]
    df.val[,topten:=ifelse(tvalue > thr, "top", "not-top")]
    df.val$title = title
    # df_subset = df.val[marker %in%  m2$marker]
    df.val[,repel:=ifelse(marker %in%  m2$marker, T, F)]
    df.val
}


plot.embedding.immunefactor.association = function(cwd, if.indexes, filename){
    
    embedding.matrix = as.matrix(tcga.phenotype[,embedding.inx,with=F])
    pheno.matrix = as.matrix(tcga.phenotype[,if.indexes,with=F])
    colnames(pheno.matrix) = gsub(colnames(pheno.matrix), pattern=".output$", replacement = "")
    # immunegene.matrix = as.matrix(tcga.phenotype[,input.genes,with=F])
    cor.curr =WGCNA::corAndPvalue(x=embedding.matrix,y=pheno.matrix)
    
    dfs.val = dfs.subset = plots.p = list()
    for (ii in seq_along(rownames(cor.curr$cor))) {
        embedding = rownames(cor.curr$cor)[ii]
        dfs.val[[embedding]] = agg.tvalue.volcano(cors=cor.curr$cor[ii,], tvalues=cor.curr$t[ii,], markers =colnames(cor.curr$cor), title=embedding)
    }
    
    dfs.val = do.call(rbind, dfs.val)
    
    
    library(dplyr)
    avi.dt =dfs.val[repel==TRUE]
    avi.dt <- avi.dt %>%
        mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40)))
    # Create wordcloud 
    library("wordcloud")
    library("RColorBrewer")
    library(ggwordcloud)
    
    set.seed(42)
    p = ggplot(
        avi.dt,
        aes(
            label = marker, size = tvalue,
            # color = factor(sample.int(10, nrow(love_words_small), replace = TRUE)), ## for random 
            angle = angle,
            color = val
            # angle = angle
        )
    ) +
        scale_color_gradient2(high="darkred", low="darkgreen", midpoint = 0) + 
        # geom_text_wordcloud_area() +
        geom_text_wordcloud_area(rm_outside = TRUE) +
        scale_size_area(max_size = 4) +
        # scale_size_area(max_size = 24) +
        theme_minimal()  +  facet_wrap(~title, ncol=8) 
    
    ggsave(p, filename = sprintf("%s/%s", cwd, filename),width = 16, height = 8)
    
}


#' Title
#'
#' @param output.dir 
#' @param sample.name 
#' @param dataset 
#' @param use.sample  NULL or vector to divide the sample. For example patient.name so that tranining and test are independent.
#' @param write.full.dataset 
#'
#' @return
#' @export
#'
#' @examples
write.dataset = function(output.dir, sample.name, dataset, use.sample=NULL, write.full.dataset=T, frac = 0.85) {
    dir.create(output.dir)
    write.table(file=paste0(output.dir, "/samples_name.txt"),x = sample.name,
                row.names = F, col.names =T,  sep="\t", quote=F )
    
    if(write.full.dataset)
        write.table(file=paste0(output.dir, "/dataset.txt"),x = dataset,
                    row.names = F, col.names =T,  sep="\t", quote=F )
    if(any(!is.null(use.sample))){
        sample.name = use.sample 
        rand_inx = sample(unique(sample.name))
        train.sample = rand_inx[1:ceiling(frac * length(rand_inx))]
        val.sample = rand_inx[ceiling(frac * length(rand_inx)+1):length(rand_inx)]
        train.inx = sample(which(sample.name %in% train.sample))
        val.inx = sample(which(sample.name %in% val.sample))
        
    }else{
        rand_inx = sample(nrow(dataset))
        train.inx = rand_inx[1:ceiling(frac * nrow(dataset))]
        val.inx = rand_inx[ceiling(frac* nrow(dataset)):nrow(dataset)]
        
    }
    
    write.table(file=paste0(output.dir, "/dataset_train.txt"),x = dataset[train.inx,],
                row.names = F, col.names =T,  sep="\t", quote=F )
    write.table(file=paste0(output.dir, "/dataset_val.txt"),x = dataset[val.inx,],
                row.names = F, col.names =T,  sep="\t", quote=F )
    
}



get_one_hot_aa = function(cdr3aa.list, width = 20, aa.uniq="ACDEFGHIKLMNPQRSTVWY"){
    xx = stringi::stri_pad(cdr3aa.list,  pad = "X", width = width, side="right")
    xx =do.call(rbind, strsplit(xx, split=""))
    xx[xx=="X"] = NA
    if(!is.null(aa.uniq)){
        aa.uniq = unlist(strsplit(aa.uniq,split = ""))
    }else{
        aa.uniq = sort(unique(c(xx)))
        aa.uniq = aa.uniq[!is.na(aa.uniq)]
        
    }
    xx = data.table(xx)
    xx =xx[, lapply(.SD, factor, levels=aa.uniq)]
    aa.one = mltools::one_hot(xx, sparsifyNAs = T)
    aa.one
}


##Check if two clusters are confounded by pre or post or patient.
text.clusters.features <- function(data, cluster,  text.cols = "condition",
                                   title="t-SNE",size=0.25,do.discrete=T, filename=NULL, normalize=TRUE, shape = 1){
    require(viridis)
    require(ggthemes)
    dt1 = as.data.frame(cluster)
    colnames(dt1) = c("V1", "V2")
    ps = list()
    for (text.col in text.cols) {
        text.col = gsub(text.col, pattern = "-", replacement = ".")
        if(normalize) data[[text.col]] = znorm(data[[text.col]])
        if(is.null(data$shape)) {
            d_cluster_1=cbind(dt1,col=data[[text.col]], shape=shape)
        }else{
            d_cluster_1=cbind(dt1,col=data[[text.col]], shape=data$shape)
        }
        d_cluster_2 = d_cluster_1[,1:2]
        title.curr = sprintf("%s_%s", title, text.col)
        p=ggplot(d_cluster_1, aes(x=V1, y=V2, colour=col), alpha=0.8, size=size) +
            geom_point(data = d_cluster_2,  colour = "grey", alpha=0.2) +
            geom_point() + 
            # geom_text(size=size,aes(label=col, colour = "grey"), alpha=0.7) +
            # geom_text( aes(label=col), alpha=0.7) +
            xlab("Dim1") + ylab("Dim2") +
            ggtitle(label = title.curr) +
            theme_light(base_size=20) +
            facet_wrap(~col)+
            theme(legend.title = element_blank())  + 
            theme(
                strip.background = element_blank(),
                strip.text.x = element_blank()
            )
        #     # theme(axis.text.x=element_blank(),
        # axis.text.y=element_blank()) 
        if(!is.null(filename)) {
            filename.curr = sprintf("%s_%s.pdf", filename, gsub(text.col, pattern="-", replacement = "_"))
            
            ggsave(file=filename.curr, p)
        }
        ps[[text.col]]  = p
        
    }
    ps
}


plot.biauc <- function(ps, dir, indexes, cell.types, mat, ncol=3) {
    require(ggrepel)
    names.ps = names(ps)
    biauc.ps = dts = list() 
    for (tt in seq_along(ps) ) {
        p = ps[[tt]]$data
        name = p$name = names.ps[tt]
        title = gsub(ps[[tt]]$label$title, pattern = " ", replacement =  ".")
        p$title = title 
        xx = title
        curr.marker = p$marker
        responder.index.curr = intersect(indexes, which(cell.types==name & response=="Responder"))
        nonresponder.index.curr = intersect(indexes, which(cell.types==name & response!="Responder"))
        auc.pos = colMeans(mat[responder.index.curr,curr.marker], na.rm = T) > colMeans(mat[nonresponder.index.curr,curr.marker], na.rm = T)
        p[,biaucs:=ifelse(auc.pos, aucs-0.5, -(aucs-0.5))]
        dts[[name]] = p 
        
        di.aucs.dt = p[label=="signature"]
        aucs.dt = p[label=="gene"]
        
        m1 = di.aucs.dt[which(aucs > 0.7)]
        if(nrow(m1) > 20) m1 = di.aucs.dt[order(aucs,decreasing=T)][1:25]
        if(nrow(m1) < 2) m1 = di.aucs.dt[order(aucs,decreasing=T)[1:5]]
        m2 = aucs.dt[which(aucs > 0.7)]
        if(nrow(m2) > 20) m2 = aucs.dt[order(aucs,decreasing=T)[1:20]]
        if(nrow(m2) < 2) m2 = aucs.dt[order(aucs,decreasing=T)[1:5]]
        
        p[,sub:=ifelse(marker %in% c(m1$marker, m2$marker), 1, 0)]
        p1 = ggplot(p, aes(x = biaucs, y = logP)) +
            geom_point(aes(color=as.factor(label), alpha = alpha)) +
            
            theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
            labs(x="AUC - 0.5", y="Significance", title=title)+
            geom_text_repel(
                data = subset(p, sub==1 ),
                aes(x = biaucs, y = logP, label = marker),
                size = 3,
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            )  
        
        biauc.ps[[name]] = p1
        filename = sprintf("%s/bi_aucs_%s.pdf",dir, title)
        ggsave(p1, file=filename, width =7, height = 7)
        
    }
    
    dts = do.call(rbind, dts)
    
    title1 = strsplit(title, split = "\\.")[[1]][1]
    auc.p = ggplot(dts, aes(x = aucs, y = logP)) +
        geom_point(aes(color=as.factor(label), alpha = alpha)) +
        
        theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
        labs(x="AUC", y="Significance", title=title1)+
        geom_text_repel(
            data = subset(dts, sub==1 ),
            aes(x = aucs, y = logP, label = marker),
            size = 3,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )   + facet_wrap(~name, ncol=ncol, scales = "free")
    filename = sprintf("%s/aucs_%s.pdf",dir, title1)
    ggsave(auc.p, file=filename, width =20, height = 15)
    
    
    biauc.p = ggplot(dts, aes(x = biaucs, y = logP)) +
        geom_point(aes(color=as.factor(label), alpha = alpha)) +
        
        theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
        labs(x="AUC - 0.5", y="Significance", title=title1)+
        geom_text_repel(
            data = subset(dts, sub==1 ),
            aes(x = biaucs, y = logP, label = marker),
            size = 3,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )   + facet_wrap(~name, ncol=ncol, scales = "free")
    filename = sprintf("%s/biaucs_%s.pdf",dir, title1)
    ggsave(biauc.p, file=filename, width =20, height = 15)
    
    list(biauc.p =biauc.p, biauc.ps=biauc.ps,auc.p=auc.p,  dts=dts)
}

get.phenotype = function(sample.subet  = NULL ) {
    dataset.sample.name = fread("/homes6/asahu/project/deeplearning/icb/data/Getz_scRNA//scrna.v4.genes/samples_name.txt")$x
    icb.phenotype = fread("/homes6/asahu/project/deeplearning/icb/data/tcga/scrna.v4.genes/tensorboardLog/nopipeline_vae_20190819-161104/epoch-142//val_prediction.csv")
    icb.phenotype = icb.phenotype[unlist(icb.phenotype$sample_name) +1]
    header = fread("/homes6/asahu/project/deeplearning/icb/data/tcga/scrna.v4.genes/tensorboardLog/nopipeline_vae_20190819-161104/best_val_0.csv",nrows = 1)
    colnames(icb.phenotype) = colnames(header)
    icb.phenotype$sample.name = dataset.sample.name
    if(!is.null(sample.subet)) icb.phenotype = icb.phenotype[match(sample.subet, sample.name)]
   icb.phenotype 
}
subset.lm.embedding= function(mat, dt, subset){
    width = dt$aa.len 
    stopifnot(sum(width)==nrow(mat))
    end = cumsum(width)
    start =  c (0, end[-length(end)]) + 1
    sub.inx = unlist(lapply(seq_along(start), function(tt) {
        out  = NULL 
        if(tt %in% subset) out = start[tt]:end[tt] 
        out
        }
        ))
    mat[sub.inx,]
}

aggregate.embedding = function(mat, dt, mode="mean"){
    width = dt$aa.len 
    stopifnot(sum(width)==nrow(mat))
    end = cumsum(width)
    start =  c (0, end[-length(end)]) + 1
    if(mode=="mean"){
        # require(IRanges)
        # by <- IRanges(start=start, width = width)
        mat.agg = do.call(rbind, lapply(seq_along(start), function(tt) colMeans(mat[start[tt]:end[tt],])))
        # mat.agg = aggregate(mat, by=by, FUN=colMeans)
    }
    if(mode=="last")
        mat.agg = mat[end,]
    if(mode=="middle"){
        middle = ceiling(start + width/2)
        mat.agg = mat[end,]
    }
    if(mode=="SIF"){
        middle = ceiling(start + width/2)
        out = calc.sif.weights(dt$aa)
        mat.agg = do.call(rbind, lapply(seq_along(start), function(tt)  out$sif.weights[start[tt]:end[tt]] %*% mat[start[tt]:end[tt],]))
        ## SVD is conducted on non-centered data 
        if(remove.pc)
          mat.agg = remove.top.pc(mat.agg, 1)
        
    }
    try(
    if(is.function(mode[[1]])){
        out = calc.positional.sif.weights(dt$aa, func = mode[[1]])
        mat.agg = do.call(rbind, lapply(seq_along(start), function(tt)  out[[tt]] %*% mat[start[tt]:end[tt],]))
        mat.agg = remove.top.pc(mat.agg, 1)
        
    })
    mat.agg 
}

remove.top.pc = function(tt, num.pc.remove=1) {
    res = prcomp(tt, rank = num.pc.remove, center = F, scale = F)
    tt - res$x %*% t(res$rotation)
}

calc.sif.weights = function(aa, a=.001){
    require(magrittr)
    xx1 = aa %>% 
        paste(collapse="") %>% 
        strsplit(split="") %>% unlist 
    xx = xx1 %>% 
        `[`(!. %in% c("", " ", ".", ",")) %>% 
        table
    xx = xx/sum(xx)
    tab = a/(a+xx)
    list(sif.weights=tab[xx1], sif.table =tab)
}

sif.smoothing = function(xx, a=.001) a/(a+xx)
    
calc.positional.sif.weights = function(aa, aa.max=26, a=.001, func=sif.smoothing){
    require(magrittr)
    xx1 = aa %>% 
        paste(collapse="") %>% 
        strsplit(split="") %>% unlist 
    xx = xx1 %>% 
        `[`(!. %in% c("", " ", ".", ",")) %>% 
        table
    all.freq = xx/sum(xx)
    all.char = names(all.freq)
    # tab.all = a/(a+xx)
    
    # forward
    xx2 = aa %>%
        strsplit(.,split="")  %>%
        sapply(.,"[", j = aa.max) 
    library(stringi)
    fwd.aa = strsplit(aa, split="")
    xx2 = stri_list2matrix( fwd.aa, byrow=TRUE)
    fwd.freq = apply(xx2, 2, function(tt) {
            uu = table(tt)
            den= max(length(tt)/3, sum(uu))
            uu = uu/den
            uu[all.char]
            }) 
    rev.aa = strsplit(stri_reverse(aa), split="")
    xx2 = stri_list2matrix(rev.aa , byrow=TRUE)
    rev.freq = apply(xx2, 2, function(tt) {
        uu = table(tt)
            den= max(length(tt)/3, sum(uu))
            uu = uu/den
            uu[all.char]
    }) 
    freq = lapply(seq_along(aa), function(tt) {
        inx.mat = cbind(match(fwd.aa[[tt]], all.char), seq_along(fwd.aa[[tt]])) 
        fwd.freq.curr = fwd.freq[inx.mat]
        inx.mat[,1] = rev(inx.mat[,1])
        rev.freq.curr = rev(rev.freq[inx.mat])
        ww = apply(cbind(fwd.freq.curr,
                   rev.freq.curr,
        all.freq[fwd.aa[[tt]]]), 1, max, na.rm=T)
        # a/(a + ww)
        func(ww)
    })
        
  freq
    
}



plot.tcr.embedding <- function(trust.curr, save.dir, title, mode="last") {
    
    # require(avinash)
    trust.curr.agg = aggregate.embedding(trust.curr, trust4.filtered, mode=mode)
    trust.curr.agg.filtered = trust.curr.agg[trust4.cdr3.filtered$currid,]

    n_neighbors =15; learning_rate =1; min_dist = .01; pca = 50
    umap.all.p = plotUMAP(data = trust.curr.agg.filtered,  
                          # umap.model = umap.all.p[[1]],
                          col=NULL, color.col=as.factor(trust4.cdr3.filtered$cdr3.type), size=1, do.discrete=T, 
                          n_neighbors = n_neighbors, learning_rate = learning_rate, min_dist = min_dist, pca=pca,
                          title= title,
                          filename=NULL , n_epochs = 100, metric = "euclidean")
    
    data.umap = data.table(umap.all.p[[1]]$embedding)
    p1 = plot.sequence.and.clustering(data.clust=data.umap, text = trust4.cdr3.filtered$aa, color.col = as.factor(trust4.cdr3.filtered$cdr3.type), num.plot.seq =100, text.size =2) 
    ggsave(filename = sprintf("%s/cdr3_seq.pdf",save.dir), p1, width = 16, height = 16)
    print(p1)
    library(RColorBrewer)
    xx1 = color.clusters.features( data=as.data.frame(phenotype_sel.curr), cluster=data.umap,  color.cols =c("assign.ident", 
                                                                                                             "assign.ident.2"),
                                   title="cell.type",size=2, filename= sprintf("%s/",save.dir), normalize=F, do.discrete=T)
    
    
    xx2 = color.clusters.features( data=as.data.frame(phenotype_sel.curr[filter1,]), cluster=data.umap[filter1,],  color.cols =c("assign.ident", 
                                                                                                                                 "assign.ident.2"),
                                   title="cell.type",size=2, filename= sprintf("%s/tcell",save.dir), normalize=F, do.discrete=T)
    print(xx2$assign.ident.2)
    
    xx3 = color.clusters.features( data=as.data.frame(icb.expression.curr[filter1,]), cluster=data.umap[filter1,],  color.cols =colnames(icb.expression.curr), title="",size=2, filename= sprintf("%s/genes/",save.dir), normalize=F, do.discrete=T)

    xx4 = color.clusters.features( data=as.data.frame(if.curr[filter1,]), cluster=data.umap[filter1,],  color.cols =colnames(if.curr)[807:1005], title="",size=2, filename= sprintf("%s/IF/",save.dir), normalize=F, do.discrete=T)
    
    list(umap.all.p = umap.all.p, xx1=xx1, xx2 = xx2, xx3 = xx3, xx4 = xx4)
}



plot.tcr.embedding.2 <- function(trust.curr, save.dir, title, mode="last", plot.seurat.marker=T) {
    dir.create(save.dir)
    
    # require(avinash)
    library(RColorBrewer)
    trust.curr = subset.lm.embedding(trust.curr, trust4.filtered, subset= trust4.cdr3.filtered.curr$currid)
    trust.curr.curr = aggregate.embedding(trust.curr, trust4.cdr3.filtered.curr, mode=mode)
 
    n_neighbors =15; learning_rate =1; min_dist = .01; pca = 50
    umap.all.p = avinash::plotUMAP(data = trust.curr.curr,  
                          col=NULL, color.col=as.factor(trust4.cdr3.filtered.curr$cdr3.type), size=1, do.discrete=T, 
                          n_neighbors = n_neighbors, learning_rate = learning_rate, min_dist = min_dist, pca=pca,
                          title= title,
                          filename=NULL , n_epochs = 100, metric = "euclidean")
    
    data.umap = data.table(umap.all.p[[1]]$embedding)
    
    p1 = plot.sequence.and.clustering(data.clust=data.umap, text = trust4.cdr3.filtered.curr$aa, color.col = as.factor(trust4.cdr3.filtered.curr$cdr3.type), num.plot.seq =100, text.size =2) 
    ggsave(filename = sprintf("%s/cdr3_seq.pdf",save.dir), p1, width = 16, height = 16)
    # print(p1)
    

    xx1 = color.clusters.features( data=as.data.frame(phenotype_sel.curr.curr), cluster=data.umap,  color.cols =c("assign.ident", "assign.ident.2"),
                                   title="cell.type",size=2, filename= sprintf("%s/",save.dir), normalize=F, do.discrete=T)
    
    uu = color.clusters.features( data=as.data.frame(trust4.cdr3.filtered.curr), cluster=data.umap,  color.cols =c("patient", "response", "treatment"),
                                   title="cell.type",size=2, filename= sprintf("%s/",save.dir), normalize=F, do.discrete=T)

    if(do.kmeans){
      
    tt = cluster.kmeans(data.umap, tsne.dir = save.dir)
    }else{
    tt = cluster.SNNGraph(data.umap, tsne.dir = save.dir)
    }
    vv = NULL
    if(plot.seurat.marker)
    vv = seurat.marker(tt$clust, tsne.dir = save.dir, sco.curr = sco.curr)
    
    # vv = seurat.marker(tt$clust[1:2268], tsne.dir = save.dir, sco.curr = sco.curr)
    # xx3 = color.clusters.features( data=as.data.frame(icb.expression.curr[filter1,]), cluster=data.umap[filter1,],  color.cols =colnames(icb.expression.curr), title="",size=2, filename= sprintf("%s/genes/",save.dir), normalize=F, do.discrete=T)
    
    # xx4 = color.clusters.features( data=as.data.frame(if.curr[filter1,]), cluster=data.umap[filter1,],  color.cols =colnames(if.curr)[807:1005], title="",size=2, filename= sprintf("%s/IF/",save.dir), normalize=F, do.discrete=T)
    
    list(umap.all.p = umap.all.p, xx1=xx1, uu = uu, tt = tt, vv = vv)
}



cluster.SNNGraph <- function(data.curr, tsne.dir="~/",  ...) {
    require(scater)
    require(scran)
    setnames(data.curr,1:2, c("UMAP1", "UMAP2"))
    data.DI = data.curr 
    
    snn.gr <- buildSNNGraph(t(data.DI), k=100, ...)
    cluster1.igraph = factor(igraph::cluster_walktrap(snn.gr, steps =4)$membership)
    data.curr$cluster = cluster1.igraph 
    p=ggplot(data.curr, aes(x=UMAP1, y=UMAP2)) +
        geom_point(size=2,aes(color=cluster), alpha=0.8) +
        guides(colour=guide_legend(override.aes=list(size=4))) +
        xlab("Dim1") + ylab("Dim2") +
        theme_light(base_size=20) +
        # scale_color_tableau() + 
        theme_classic() +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank()) 
    # print(p)
    filename = sprintf("%s/SNN_clustering.pdf", tsne.dir) 
    ggsave(file=filename, p)
    list(clust=data.curr$cluster, p=p)
}

cluster.kmeans <- function(data.curr, centers=3, tsne.dir="~/", ...) {
  setnames(data.curr,1:2, c("UMAP1", "UMAP2"))
  set.seed(100)
  
  embedding.mat  = data.curr[,.SD, .SDcols=names(data.curr) %like% "^UMAP"]
  data.DI = embedding.mat 
  clust.kmeans <- kmeans(data.DI, centers=centers, nstart = 25)
  # table(clust.kmeans$cluster)
  
  data.curr$cluster1 = (clust.kmeans$cluster)
  data.curr[,cluster:= as.factor(cluster1)]
  p=ggplot(data.curr, aes(x=UMAP1, y=UMAP2)) +
    geom_point(size=2,aes(color=cluster), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=4))) +
    xlab("Dim1") + ylab("Dim2") +
    # ggtitle(label = title) +
    theme_light(base_size=20) +
    scale_color_tableau() + 
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
  
    filename = sprintf("%s/kmean_clustering.pdf", tsne.dir) 
    ggsave(file=filename, p)

  list(clust=data.curr$cluster, p=p)
}
seurat.marker <- function(clust, tsne.dir, sco.curr){
    require(Seurat)
    require(pROC)
    
    sco.curr@meta.data$deepImmune.clust = clust
    Idents(object = sco.curr) <- sco.curr@meta.data$deepImmune.clust
    # sco.curr.markers <- FindAllMarkers(sco.curr, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.3)
    sco.curr.markers <- FindAllMarkers(sco.curr)
    # sco.markers = find.markers(sco.curr,7,thresh.use = 2,test.use = "roc")
    # top2 = sco.curr.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    # top10 <- sco.curr.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% subset(p_val < 1E-3)
    top10 <- sco.curr.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)%>% subset(p_val < 1E-3)
    vlp.p = VlnPlot(sco.curr,features=top10$gene,  slot = "counts", log = TRUE)
    ggsave(file=sprintf("%s/vlnplot_markers.pdf", tsne.dir), vlp.p, width=20, height=20)
    ## heatmap plot 
    p = DoHeatmap(sco.curr, features = top10$gene) + NoLegend()
    ggsave(file=sprintf("%s/heatmap_markers.pdf", tsne.dir), p, width=16, height=16)
    
    top10 <- sco.curr.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)%>% subset(p_val < 1E-3)
    top10 = data.table(top10)
    write.table(file=sprintf("%s/top_markers.txt", tsne.dir), x=top10, row.names = F, col.names = T, quote = F)
    top10 = top10[gene %in% genentech.feat$symbol]
    top10$cluster = as.numeric(top10$cluster)
    markers.auc = list()
    for(ii in unique(top10$cluster)){
        genes.curr = top10[cluster==ii]$gene
        out = lapply(genes.curr, function(gene.curr) {
            exp.curr =genentech.exp[match(gene.curr, genentech.feat$symbol),]
            calc.stat.new(genentech.response, exp.curr)
        })
        out.dt = data.table(do.call(rbind, out))
        out.dt$genes = genes.curr
        out.dt$clust = ii
        if(nrow(out.dt) >= 1) 
            markers.auc= append(markers.auc, list(out.dt))
    }
    
    markers.auc.dt = do.call(rbind, markers.auc)
    markers.auc.dt = markers.auc.dt[order(V1,decreasing =T)]
    # proabably use SSGSEA to estiamte the fraction of tumors 
    require(ggrepel)
    p = ggplot(data=markers.auc.dt, aes(x=as.factor(clust), y = V1)) +
        geom_boxplot() + geom_point() + 
        geom_boxplot(outlier.shape=10, outlier.size=8)  +
        geom_point(position="jitter", size=2) + 
        geom_text_repel(
            data = markers.auc.dt[V1>0.6],
            aes(x=as.factor(clust),  y = V1, label = genes),
            size = 3,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )  
    ggsave(file=sprintf("%s/marker_genentech_response.pdf", tsne.dir), p, width=16, height=10)
    NULL
}




translate.seq = function(tt) {
  ifelse(stringr::str_length(tt) %% 3 ==0, 
         as.character(Biostrings::translate(Biostrings::DNAStringSet(tt),if.fuzzy.codon = "X")), NA)
}
get.aminoacid = function(CDR3, V , J){
  missing = "X"
  out = NA
  str.len = stringr::str_length(CDR3)
  if(V!="*" & J=="*"){
    CDR3 = substr(CDR3,1, floor(str.len/3) *3)
    out = translate.seq(CDR3)
    out = stringr::str_c( out, missing)
  }else if(V=="*" & J!="*"){
    # CDR3 = stringi::stri_reverse(CDR3)
    CDR3 = substr(CDR3,1 + ( str.len %%3), str.len)
    out = translate.seq(CDR3)
    out = stringr::str_c(missing, out)
    # out = stringi::stri_reverse(out)
  }else if(V!="*" & J!="*" & str.len %%3 == 0){
    out = translate.seq(CDR3)
  }
  out
}
get.aminoacid.vec = Vectorize(get.aminoacid)

my.translate = function(dt, name=NULL, min.len= 6, max.len = 25){
  out = NULL 
  if(nrow(dt)){
    cols = c("consesus_id", "the_index_of_CDR3_in_that_consensus", "V", "J", "C", "CDR1", "CDR2", "CDR3", "score_of_CDR3", "abundance")
    setnames(dt,1:10, cols)
    dt$name = name 
    dt[,aa:=get.aminoacid.vec(CDR3, V , J)]
    dt[,aa.len:=stringr::str_length(aa)]
    dt[,is.alpha:=ifelse(grepl("TRA", V)|grepl("TRA", J)| grepl("TRA", C), T, F)]
    dt[,abs.score:= abundance * (score_of_CDR3+.001)]
    if(!is.na(min.len)) dt = dt[aa.len>= min.len]
    if(!is.na(max.len)) dt = dt[aa.len< max.len]
    dt.alpha = dt[is.alpha==T][abs.score == max(abs.score)]
    if(nrow(dt.alpha) > 1) dt.alpha = dt.alpha[which.min(aa.len)]
    dt.beta = dt[is.alpha==F][abs.score == max(abs.score, na.rm = T)]
    if(nrow(dt.beta) > 1) dt.beta = dt.beta[which.min(aa.len)]
    out = rbind(dt.alpha, dt.beta)
  }
  out
}



get.trust4.cdr3 <- function() {
  


library(stringr)
library(stringi)
load("/liulab/asahu/data/ssgsea/xiaoman/getz/trust4.out.V2.RData")
names(trust4.out) = c("cdr3", "annot")
substr.id = function(tt){
  paste(strsplit(tt, split="_")[[1]][1:3], collapse = "_")
}

trust4.cellid = sapply(names(trust4.out[[1]]), substr.id)
tcell.id = sapply(rownames(icb.expression.matched),  substr.id)
sum(trust4.cellid %in% tcell.id) ## only 3196 
only.tcells = setdiff(tcell.id, trust4.cellid)
common.cellids = intersect(trust4.cellid,tcell.id)
cdr3.tcell = trust4.out$cdr3[match(common.cellids, trust4.cellid)]
sample.name.map =data.table(sample.name.short=common.cellids)
sample.name.map$sample.name = rownames(icb.expression.matched)[match(common.cellids, tcell.id)]
sample.name.map$trust.id = names(trust4.out[[1]])[match(common.cellids, trust4.cellid)]

trust4.cdr3  =do.call(rbind,  lapply(seq_along(trust4.out$cdr3), function(tt) {
  dt = trust4.out$cdr3[[tt]]
  if(nrow(dt) ==0 ){ 
    dt = NULL
  }else{
    dt$cellid = names(trust4.out$cdr3)[[tt]]
  }
  dt
}))
min.len = 7*3
max.len =25*3
cols = c("consesus_id", "the_index_of_CDR3_in_that_consensus", "V", "J", "C", "CDR1", "CDR2", "CDR3", "score_of_CDR3", "abundance")
setnames(trust4.cdr3,1:10, cols)
trust4.cdr3$sample.name.short = sample.name.map[match(trust4.cdr3$cellid,trust.id)]$sample.name.short
trust4.cdr3$sample.name = sample.name.map[match(trust4.cdr3$cellid,trust.id)]$sample.name
trust4.cdr3[,str.len:=stringr::str_length(CDR3)]
trust4.cdr3[,is.alpha:=ifelse(grepl("TRA", V)|grepl("TRA", J)| grepl("TRA", C), T, F)]
trust4.cdr3[,is.beta:=ifelse(grepl("TRB", V)|grepl("TRB", J)| grepl("TRB", C), T, F)]
trust4.cdr3[,abs.score:= abundance * (score_of_CDR3+.001)]
# trust4.cdr3[, cellid.short:=sapply(cellid, substr.id)]
if(!is.na(min.len)) trust4.cdr3 = trust4.cdr3[str.len>= min.len]
if(!is.na(max.len)) trust4.cdr3 = trust4.cdr3[str.len<= max.len]
trust4.cdr3
}