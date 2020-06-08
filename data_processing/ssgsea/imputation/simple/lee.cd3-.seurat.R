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
lee.sco <-  as.Seurat(sce[,sce$CD3.status=="CD3-"], data=NULL, project="LEE") %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE)
lee.sco = lee.sco %>% 
    RunPCA(pc.genes = lee.sco@var.genes, npcs = 20, verbose = FALSE)

rm(sce); gc()



options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = lee.sco, reduction = "pca", pt.size = .1, group.by = "patient.name", do.return = TRUE)
# p2 <- VlnPlot(object = lee.sco, features = "PC_1", group.by = "patient.name", do.return = TRUE, pt.size = .1)
pdf(".figs/lee.cd3neg/lee.sco.without.harmony.pdf")
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

pdf(".figs/lee.cd3neg/lee.sco.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

lee.sco[["umapunorm"]] = lee.sco@reductions$pca@cell.embeddings %>% 
umap(., n_threads=32, metric="cosine", pca=NULL, n_neighbors=50) %>%
set_rownames(colnames(lee.sco)) %>%
set_colnames(paste0("umapunorm_", 1:2)) %>%
CreateDimReducObject(embeddings = ., key = "umapunorm_", assay = DefaultAssay(lee.sco))

DimPlot(lee.sco, reduction = "umapunorm", group.by = "patient.name", pt.size = 0.5)%>%
ggsave(filename=".figs/lee.cd3neg/lee.sco.umap.without.harmony.pdf", plot=.)

lee.sco <- lee.sco %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

DimPlot(lee.sco, reduction = "umap", group.by = "patient.name", pt.size = .1, split.by = 'patient') %>%
ggsave(".figs/lee.cd3neg/lee.sco.harmony.umap.patient.pdf", plot=.)

DimPlot(lee.sco, reduction = "umap", group.by = "patient.name", label = TRUE, pt.size = .1) %>%
ggsave(".figs/lee.cd3neg/lee.sco.harmony.umap.pdf", plot=.)

DimPlot(lee.sco, reduction = "umap", label = TRUE, pt.size = .1) %>%
ggsave(".figs/lee.cd3neg/lee.sco.harmony.umap.cluster.pdf", plot=.)

## 
p1 = DimPlot(lee.sco, reduction = "umap", group.by = "response", pt.size = .1, split.by = 'response')
pdf(".figs/lee.cd3neg/lee.sco.harmony.umap.response.pdf")
plot_grid(p1)
dev.off()
saveRDS(file="~/liulab_home/data/single_cell/scrna.dataset.deg.RData",scrna.dataset.list)