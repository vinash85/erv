## convert to Seurat object and preprocess
library(data.table)
library(magrittr)
library(ggplo2)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(cowplot)
library(harmony)
library(uwot)
library(parallel)
library(EnhancedVolcano)
library(tidyverse)
setwd("~/project/deeplearning/icb/deepImmune/data_processing/ssgsea/imputation/simple")
sce = readRDS("~/liulab_home/data/single_cell/Lee_data/lee.sce.RDS")
lee.sco <-  as.Seurat(sce[,sce$CD3.status=="CD3+"], data=NULL, project="LEE") %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE)
lee.sco = lee.sco %>% 
    RunPCA(pc.genes = lee.sco@var.genes, npcs = 20, verbose = FALSE)

rm(sce); gc()



options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lee.sco, reduction = "pca", pt.size = .1, group.by = "patient.name", do.return = TRUE)
# p2 <- VlnPlot(object = lee.sco, features = "PC_1", group.by = "patient.name", do.return = TRUE, pt.size = .1)
pdf(".figs/lee/lee.sco.without.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

lee.sco <- lee.sco %>% 
    RunHarmony("patient.name")


harmony_embeddings <- Embeddings(lee.sco, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lee.sco, reduction = "harmony", pt.size = .1, group.by = "patient.name", do.return = TRUE)
# p2 <- VlnPlot(object = lee.sco, features = "harmony_1", group.by = "patient.name", do.return = TRUE, pt.size = .1)

pdf(".figs/lee/lee.sco.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

lee.sco[["umapunorm"]] = lee.sco@reductions$pca@cell.embeddings %>% 
umap(., n_threads=32, metric="cosine", pca=NULL, n_neighbors=50) %>%
set_rownames(colnames(lee.sco)) %>%
set_colnames(paste0("umapunorm_", 1:2)) %>%
CreateDimReducObject(embeddings = ., key = "umapunorm_", assay = DefaultAssay(lee.sco))

DimPlot(lee.sco, reduction = "umapunorm", group.by = "patient.name", pt.size = 0.5)%>%
ggsave(filename=".figs/lee/lee.sco.umap.without.harmony.pdf", plot=.)

lee.sco <- lee.sco %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

DimPlot(lee.sco, reduction = "umap", group.by = "patient.name", pt.size = .1, split.by = 'patient') %>%
ggsave(".figs/lee/lee.sco.harmony.umap.patient.pdf", plot=.)

DimPlot(lee.sco, reduction = "umap", group.by = "patient.name", label = TRUE, pt.size = .1) %>%
ggsave(".figs/lee/lee.sco.harmony.umap.pdf", plot=.)

DimPlot(lee.sco, reduction = "umap", label = TRUE, pt.size = .1) %>%
ggsave(".figs/lee/lee.sco.harmony.umap.cluster.pdf", plot=.)

## 
p1 = DimPlot(lee.sco, reduction = "umap", group.by = "response", pt.size = .1, split.by = 'response')
pdf(".figs/lee/lee.sco.harmony.umap.response.pdf")
plot_grid(p1)
dev.off()

saveRDS(file="~/liulab_home/data/single_cell/Lee_data/lee.seurat.cd3.RDS", object=lee.sco)

## For each cluster find differential expression 
# Define models
library(lme4)
library(tidyverse)
f1 <- 'expr ~ offset(log_total_count) + (1 | patient)'
f2 <- 'expr ~ offset(log_total_count) + (1 | binaryResponse/patient)'
eval.lmer.poisson = function(tmp_df, val){
	tryCatch({
		if(sum(val>0) <=1 ) stop()
		tmp_df = tmp_df %>%
		mutate(expr=val)
		mres1 <- glmer(f1, data = tmp_df, family = poisson())
		mres2 <- glmer(f2, data = tmp_df, family = poisson())
		anova_res <- anova(mres1, mres2)
		pval <- anova_res$`Pr(>Chisq)`[2]
		effect_size <- ranef(mres2)$binaryResponse[,1][1] - ranef(mres2)$binaryResponse[,1][2]

		c(pval = pval, effect_size = effect_size)
	}, error = function(e) c(pval = NA, effect_size = NA))
}

# for each cluster 
cluster.sel = names(which(table(lee.sco$seurat_clusters)>1000))
# all.clusters= unique(lee.sco$seurat_clusters)
pre.diff.dt =list()
for (clust in cluster.sel) {
	lee.sco.pre =lee.sco[, (lee.sco$Treatment.Cycle=="C1D1") & (lee.sco$seurat_clusters==clust) & (lee.sco$response %in% c("CR", "PR", "PD"))]
	exp.curr = lee.sco.pre@assays$RNA@counts
	exp.curr = exp.curr[rowSums(exp.curr)>0, ]
	meta.dt = lee.sco.pre@meta.data %>%
	as.data.table() %>%
	.[,.(binaryResponse=ifelse(response %in% c("CR", "PR"),1 ,0) , patient=patient.name)] %>%
	.[,total_count:=colSums(exp.curr, na.rm=T)] %>%
	.[,log_total_count:=log(total_count)] 

	out = mclapply(seq(nrow(exp.curr)), function(tt) 
		eval.lmer.poisson(meta.dt, exp.curr[tt,]),
		mc.cores=48
		)
# out = mclapply(seq(100), function(tt) 
#     eval.lmer.poisson(meta.dt, exp.curr[tt,]),
#     mc.cores=32
    # )
	clust.label = sprintf("clust%s", clust) 
	pre.diff.dt[[clust]] = do.call(rbind, out)%>% 
	as.data.table() %>%
	.[,gene:=rownames(exp.curr)] %>%
	.[!is.na(pval)] 


	EnhancedVolcano(pre.diff.dt[[clust]],
		lab = pre.diff.dt[[clust]]$gene,
		x = 'effect_size',
		y = 'pval',
		pCutoff = 1e-4,
		FCcutoff = 1,
                             # ylim = c(0,5.2),
                             # xlim = c(-1, 1),
                             # pointSize = 4.0,
		pointSize = c(ifelse(pre.diff.dt[[clust]]$pval< 1E-4, 1, 0.2)),
		labSize = 4.0,
		legend=c('NS','Log (base 2) fold-change','Adj.P value',
			'Adj.P value & Log (base 2) fold-change'),
		legendPosition = 'right',
		legendLabSize = 8,
		legendIconSize = 4.0,
		drawConnectors = TRUE,
		widthConnectors = 0.2,
		colAlpha = 0.8,
		colConnectors = 'grey30'
		) %>%
	ggsave(sprintf(".figs/lee/lee.%s.volcano.pre.pdf", clust), .)

}

save(file=".figs/lee/lee.cd3.pre.deg.RData",pre.diff.dt)


post.diff.dt =list()
for (clust in cluster.sel) {
	lee.sco.post =lee.sco[, (lee.sco$Treatment.Cycle=="C4D1") & (lee.sco$seurat_clusters==clust) & (lee.sco$response %in% c("CR", "PR", "PD"))]
	exp.curr = lee.sco.post@assays$RNA@counts
	exp.curr = exp.curr[rowSums(exp.curr)>0, ]
	meta.dt = lee.sco.post@meta.data %>%
	as.data.table() %>%
	.[,.(binaryResponse=ifelse(response %in% c("CR", "PR"),1 ,0) , patient=patient.name)] %>%
	.[,total_count:=colSums(exp.curr, na.rm=T)] %>%
	.[,log_total_count:=log(total_count)] 

	out = mclapply(seq(nrow(exp.curr)), function(tt) 
		eval.lmer.poisson(meta.dt, exp.curr[tt,]),
		mc.cores=48
		)
# out = mclapply(seq(100), function(tt) 
#     eval.lmer.poisson(meta.dt, exp.curr[tt,]),
#     mc.cores=32
    # )
	clust.label = sprintf("clust%s", clust) 
	post.diff.dt[[clust]] = do.call(rbind, out)%>% 
	as.data.table() %>%
	.[,gene:=rownames(exp.curr)] %>%
	.[!is.na(pval)] 


	EnhancedVolcano(post.diff.dt[[clust]],
		lab = post.diff.dt[[clust]]$gene,
		x = 'effect_size',
		y = 'pval',
		pCutoff = 1e-4,
		FCcutoff = 1,
                             # ylim = c(0,5.2),
                             # xlim = c(-1, 1),
                             # pointSize = 4.0,
		pointSize = c(ifelse(post.diff.dt[[clust]]$pval< 1E-4, 1, 0.2)),
		labSize = 4.0,
		legend=c('NS','Log (base 2) fold-change','Adj.P value',
			'Adj.P value & Log (base 2) fold-change'),
		legendPosition = 'right',
		legendLabSize = 8,
		legendIconSize = 4.0,
		drawConnectors = TRUE,
		widthConnectors = 0.2,
		colAlpha = 0.8,
		colConnectors = 'grey30'
		) %>%
	ggsave(sprintf(".figs/lee/lee.%s.volcano.post.pdf", clust), .)

}

save(file=".figs/lee/lee.cd3.post.deg.RData",post.diff.dt)




## global differences pretreatment 

lee.sco.pre =lee.sco[, (lee.sco$Treatment.Cycle=="C1D1") &  & (lee.sco$response %in% c("CR",  "PD"))]
exp.curr = lee.sco.pre@assays$RNA@counts
exp.curr = exp.curr[rowSums(exp.curr)>0, ]
meta.dt = lee.sco.pre@meta.data %>%
as.data.table() %>%
.[,.(binaryResponse=ifelse(response %in% c("CR", "PR"),1 ,0) , patient=patient.name)] %>%
.[,total_count:=colSums(exp.curr, na.rm=T)] %>%
.[,log_total_count:=log(total_count)] 

out = mclapply(seq(nrow(exp.curr)), function(tt) 
	eval.lmer.poisson(meta.dt, exp.curr[tt,]),
	mc.cores=48
	)
global.pre = do.call(rbind, out)%>% 
as.data.table() %>%
.[,gene:=rownames(exp.curr)] %>%
.[!is.na(pval)] 


EnhancedVolcano(global.pre,
	lab = global.pre$gene,
	x = 'effect_size',
	y = 'pval',
	pCutoff = 1e-4,
	FCcutoff = 1,
                         # ylim = c(0,5.2),
                         # xlim = c(-1, 1),
                         # pointSize = 4.0,
	pointSize = c(ifelse(global.pre$pval< 1E-4, 1, 0.2)),
	labSize = 4.0,
	legend=c('NS','Log (base 2) fold-change','Adj.P value',
		'Adj.P value & Log (base 2) fold-change'),
	legendPosition = 'right',
	legendLabSize = 8,
	legendIconSize = 4.0,
	drawConnectors = TRUE,
	widthConnectors = 0.2,
	colAlpha = 0.8,
	colConnectors = 'grey30'
	) %>%
ggsave(".figs/lee/lee.global.volcano.pre.offset.pdf", .)

save(file=".figs/lee/lee.cd3.global.pre.deg.offset.RData",global.pre)




## global differences posttreatment 

lee.sco.post =lee.sco[, (lee.sco$Treatment.Cycle=="C4D1") &  & (lee.sco$response %in% c("CR",  "PD"))]
exp.curr = lee.sco.post@assays$RNA@counts
exp.curr = exp.curr[rowSums(exp.curr)>0, ]
meta.dt = lee.sco.post@meta.data %>%
as.data.table() %>%
.[,.(binaryResponse=ifelse(response %in% c("CR", "PR"),1 ,0) , patient=patient.name)] %>%
.[,total_count:=colSums(exp.curr, na.rm=T)] %>%
.[,log_total_count:=log(total_count)] 

out = mclapply(seq(nrow(exp.curr)), function(tt) 
	eval.lmer.poisson(meta.dt, exp.curr[tt,]),
	mc.cores=48
	)
global.post = do.call(rbind, out)%>% 
as.data.table() %>%
.[,gene:=rownames(exp.curr)] %>%
.[!is.na(pval)] 


EnhancedVolcano(global.post,
	lab = global.post$gene,
	x = 'effect_size',
	y = 'pval',
	pCutoff = 1e-4,
	FCcutoff = 1,
                         # ylim = c(0,5.2),
                         # xlim = c(-1, 1),
                         # pointSize = 4.0,
	pointSize = c(ifelse(global.post$pval< 1E-4, 1, 0.2)),
	labSize = 4.0,
	legend=c('NS','Log (base 2) fold-change','Adj.P value',
		'Adj.P value & Log (base 2) fold-change'),
	legendPosition = 'right',
	legendLabSize = 8,
	legendIconSize = 4.0,
	drawConnectors = TRUE,
	widthConnectors = 0.2,
	colAlpha = 0.8,
	colConnectors = 'grey30'
	) %>%
ggsave(".figs/lee/lee.global.volcano.post.offset.pdf", .)

save(file=".figs/lee/lee.cd3.global.post.deg.offset.RData",global.post)