
library(Matrix)
library(Seurat)
library(data.table)
library(magrittr)
library(cowplot)
library(harmony)
library(uwot)
library(parallel)
exp.file = list.files(path = "~/liulab_home/data/single_cell/GSE145281/", pattern="*raw.txt.gz", full.names=T) %>%
grep(pattern="R[1-5]", value=T)
GSE145281.list = lapply(exp.file, function(tt){
    aa = fread(tt)
    aa[,-1,with=F] %>% 
    as.matrix() %>%
    t()%>%
    set_colnames(aa[[1]])%>%
    Matrix(data=., sparse=T)%>%
    list(exp=.,samples=rownames(.), genes=colnames(.))
}
    )
GSE145281.exp.mat= sapply(GSE145281.list, function(tt) tt$exp, simplify=F) %>%
do.call(rbind,.)
GSE145281.meta = lapply(seq(length(exp.file)),  function(tt) {
    patient.curr= basename(exp.file[[tt]]) %>%
    strsplit(x=., split="_")%>% 
    unlist() %>%.[2]
    response = substring(patient.curr,1,1) %>%
    equals("R")+0
    data.table(samples=GSE145281.list[[tt]]$samples, patient=patient.curr, response=response)
}) %>% do.call(rbind, .)

GSE145281.umap.model.exp <- umap(as.matrix(GSE145281.exp.mat), pca=50, n_threads=32)
d_GSE145281.umap = cbind(GSE145281.meta[,.(patient=patient, response=response)], data.table(GSE145281.umap.model.exp)) %>%
.[,label:=as.factor(patient)]
levels(d_GSE145281.umap$label) = LETTERS[seq(length(unique(d_GSE145281.umap$label)))]
p = ggplot(d_GSE145281.umap, aes(x = V1, y = V2)) + 
  geom_point(aes(color = as.factor(label)), alpha = 0.5, size=0.5) +
  # geom_text(aes(color = as.factor(response), label = label), alpha = 0.7, size=2) +
  guides(colour = guide_legend(override.aes = list(size = 1))) + 
  xlab("umap_1") + ylab("umap_2")  +
  theme_classic()
ggsave(".figs/GSE145281.umap.expression.pdf", p)


# patient specific correction 
gse145281 <- CreateSeuratObject(counts = t(GSE145281.exp.mat), project = "GSE145281", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)

GSE145281.meta.match = GSE145281.meta[match(rownames(gse145281@meta.data), samples)]
gse145281@meta.data %<>% cbind(., GSE145281.meta.match[,.(patient, response)]) %>%
as.data.frame()  

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gse145281, reduction = "pca", pt.size = .1, group.by = "patient", do.return = TRUE)
p2 <- VlnPlot(object = gse145281, features = "PC_1", group.by = "patient", do.return = TRUE, pt.size = .1)
pdf(".figs/gse145281.without.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

gse145281 <- gse145281 %>% 
    RunHarmony("patient")


harmony_embeddings <- Embeddings(gse145281, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gse145281, reduction = "harmony", pt.size = .1, group.by = "patient", do.return = TRUE)
p2 <- VlnPlot(object = gse145281, features = "harmony_1", group.by = "patient", do.return = TRUE, pt.size = .1)

pdf(".figs/gse145281.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

gse145281[["umapunorm"]] = gse145281@reductions$pca@cell.embeddings %>% 
umap(., n_threads=32, metric="cosine", pca=NULL, n_neighbors=50) %>%
set_rownames(colnames(gse145281)) %>%
set_colnames(paste0("umapunorm_", 1:2)) %>%
CreateDimReducObject(embeddings = ., key = "umapunorm_", assay = DefaultAssay(gse145281))

DimPlot(gse145281, reduction = "umapunorm", group.by = "patient", pt.size = 0.5)%>%
ggsave(filename=".figs/gse145281.umap.without.harmony.pdf", plot=.)

gse145281 <- gse145281 %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

DimPlot(gse145281, reduction = "umap", group.by = "patient", pt.size = .1, split.by = 'patient') %>%
ggsave(".figs/gse145281.harmony.umap.patient.pdf", plot=.)

DimPlot(gse145281, reduction = "umap", group.by = "patient", label = TRUE, pt.size = .1) %>%
ggsave(".figs/gse145281.harmony.umap.pdf", plot=.)

DimPlot(gse145281, reduction = "umap", label = TRUE, pt.size = .1) %>%
ggsave(".figs/gse145281.harmony.umap.cluster.pdf", plot=.)

## 
p1 = DimPlot(gse145281, reduction = "umap", group.by = "response", pt.size = .1, split.by = 'response')
pdf(".figs/gse145281.harmony.umap.response.pdf")
plot_grid(p1)
dev.off()

## Implementation of psuedo-bulk method from : https://www.biorxiv.org/content/10.1101/2020.04.01.019851v1 

```{r}
gse145281.clust = gse145281[,gse145281$seurat_clusters==2]
patient.cell.one.hot= gse145281.clust@meta.data %>% 
data.table(ID=rownames(.)) %>%
.[,.(ID, patient=as.factor(patient))] %>%
.[,.(ID, patient)]  %>%
mltools::one_hot(.) %>%
.[,-1,with=F] %>%
as.matrix() 

psuedo.bulk = gse145281.clust@assays$RNA@counts %*% patient.cell.one.hot

cells.1 = grep("_R", colnames(psuedo.bulk), value=T)
cells.2 = grep("_NR", colnames(psuedo.bulk), value=T)
deseq.out = DESeq2DETest(data.use=psuedo.bulk, cells.1=cells.2, cells.2=cells.1)
deseq.dt = deseq.out %>%
as.data.frame() %>%
mutate(gene=rownames(.)) %>%
data.table() 
deseq.dt = deseq.dt[order(padj)]
deseq.dt.matched = deseq.dt[match(markers.merge$gene,gene)]

markers.merge$pseudoBulk=deseq.dt[match(markers.merge$genes,gene)]$stat
cor.test(out.final.clust2$t, deseq.dt[match(out.final.clust2$gene,gene)]$stat, method="spearman")
cor(markers.merge[,12:17,with=F], use = "pairwise.complete.obs")
## check top 25 candidates of each 
# effect size and p-values
for (ii in seq(12,17)) {
    label = colnames(markers.merge)[ii]
    genes.curr = order(abs(markers.merge[[label]]), decreasing=T)[1:16] %>%
    markers.merge[.] %>%
    .[["genes"]]
    VlnPlot(gse145281.clust, genes.curr, group.by = 'patient') %>%
    ggsave(file=sprintf(".figs/gse145281.clust2.%s.pdf", label), ., width=15, height=10)
}

genes.curr = deseq.dt[order(abs(stat), decreasing=T)[1:16]]$gene

VlnPlot(gse145281.clust, genes.curr, group.by = 'patient') %>%
ggsave(file=sprintf(".figs/gse145281.clust2.%s.nonvariablegenes.pdf", "psuedoBulk"), ., width=15, height=10)



## singleR annotation 
sco = readRDS(file="~/liulab_home/data/single_cell/GSE145281/seurat.RDS") 
singleR.annotations = list()
library(SingleR)
hpca.se <- MonacoImmuneData()
singleR.annotations$MonacoImmuneData.fine <- SingleR(test = sco@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.fine)
singleR.annotations$MonacoImmuneData.main <- SingleR(test = sco@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main)

hpca.se <- HumanPrimaryCellAtlasData()
singleR.annotations$HumanPrimaryCellAtlasData.fine <- SingleR(test = sco@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.fine)
singleR.annotations$HumanPrimaryCellAtlasData.main <- SingleR(test = sco@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main)
saveRDS(file="~/liulab_home/data/single_cell/GSE145281/singleR.annotations.RDS", singleR.annotations) 

gse145281$HumanPrimaryCellAtlasData.fine=singleR.annot$HumanPrimaryCellAtlasData.fine[colnames(gse145281),]$pruned.labels
gse145281$HumanPrimaryCellAtlasData.main=singleR.annot$HumanPrimaryCellAtlasData.main[colnames(gse145281),]$pruned.labels
gse145281$MonacoImmuneData.fine=singleR.annot$MonacoImmuneData.fine[colnames(gse145281),]$pruned.labels
gse145281$MonacoImmuneData.main=singleR.annot$MonacoImmuneData.main[colnames(gse145281),]$pruned.labels