library(Seurat)
library(ggplot2)

ESRRA_all <- as.character(read.csv('ESRRA_taget_filtered.csv')[,1])
ESRRA_oxphos <- as.character(read.csv('ESRRA_oxphos_gene_with_highRP.csv')[,1])

BCC <- readRDS('/homes/cwang/projects/MAESTRO/analysis/ALIGN/MAESTRO/GSE129785_BCC_All_MAESTRO_SeuratObj.rds')
BCC_RNA <- subset(BCC, cells = rownames(BCC@meta.data[which(BCC@meta.data[,'tech']=='RNA'),]))
BCC_ATAC <- subset(BCC, cells = rownames(BCC@meta.data[which(BCC@meta.data[,'tech']=='ATAC'),]))

#' Average the expression level/regulatory potential score of ESRRA_all targets and ESRRA_oxphos targets.
BCC_RNA@meta.data$ESRRA_all <- colMeans(x = as.matrix(GetAssayData(BCC_RNA))[intersect(ESRRA_all,rownames(BCC_RNA)), ], na.rm = TRUE)
BCC_RNA@meta.data$ESRRA_oxphos <- colMeans(x = as.matrix(GetAssayData(BCC_RNA))[intersect(ESRRA_oxphos,rownames(BCC_RNA)), ], na.rm = TRUE)
BCC_ATAC@meta.data$ESRRA_all <- colMeans(x = as.matrix(GetAssayData(BCC_ATAC))[intersect(ESRRA_all,rownames(BCC_ATAC)), ], na.rm = TRUE)
BCC_ATAC@meta.data$ESRRA_oxphos <- colMeans(x = as.matrix(GetAssayData(BCC_ATAC))[intersect(ESRRA_oxphos,rownames(BCC_ATAC)), ], na.rm = TRUE)

p <- VlnPlot(BCC_RNA, assay="RNA", features = "ESRRA",  pt.size = 0, group.by = "assign.ident") + labs(title = "Expression of ESRRA in BCC",y = "Expression level (logTPM/10)") + NoLegend()
ggsave(file.path(paste0(BCC_RNA@project.name, "_ESRRA_RNA.pdf")), p, width=8, height=5, useDingbats=FALSE)   
p <- VlnPlot(BCC_ATAC, assay="ACTIVITY", features = "ESRRA",  pt.size = 0, group.by = "assign.ident") + labs(title = "Regulation of ESRRA in BCC",y = "Regulatory potential (log2)") + NoLegend()
ggsave(file.path(paste0(BCC_ATAC@project.name, "_ESRRA_ATAC.pdf")), p, width=8, height=5, useDingbats=FALSE) 
p <- VlnPlot(BCC_RNA, assay="RNA", features = "ESRRA_all",  pt.size = 0, group.by = "assign.ident") + labs(title = "Expression of ESRRA all targets in BCC",y = "Expression level (logTPM/10)") + NoLegend()
ggsave(file.path(paste0(BCC_RNA@project.name, "_ESRRA_all_RNA.pdf")), p, width=8, height=5, useDingbats=FALSE)   
p <- VlnPlot(BCC_ATAC, assay="ACTIVITY", features = "ESRRA_all",  pt.size = 0, group.by = "assign.ident") + labs(title = "Regulation of ESRRA all targets in BCC",y = "Regulatory potential (log2)") + NoLegend()
ggsave(file.path(paste0(BCC_ATAC@project.name, "_ESRRA_all_ATAC.pdf")), p, width=8, height=5, useDingbats=FALSE) 
p <- VlnPlot(BCC_RNA, assay="RNA", features = "ESRRA_oxphos",  pt.size = 0, group.by = "assign.ident") + labs(title = "Expression of ESRRA oxphos targets in BCC",y = "Expression level (logTPM/10)") + NoLegend()
ggsave(file.path(paste0(BCC_RNA@project.name, "_ESRRA_oxphos_RNA.pdf")), p, width=8, height=5, useDingbats=FALSE)   
p <- VlnPlot(BCC_ATAC, assay="ACTIVITY", features = "ESRRA_oxphos",  pt.size = 0, group.by = "assign.ident") + labs(title = "Regulation of ESRRA oxphos targets in BCC",y = "Regulatory potential (log2)") + NoLegend()
ggsave(file.path(paste0(BCC_ATAC@project.name, "_ESRRA_oxphos_ATAC.pdf")), p, width=8, height=5, useDingbats=FALSE)

#' Scanning the enriched regulators in scATAC-seq clusters using GIGGLE.
BCC_GIGGLE <- read.table('/homes/cwang/projects/MAESTRO/analysis/ATAC/cluster/MAESTRO/GSE129785_BCC_All_TF.GIGGLE_ChIPseq.txt',check.names=F)
BCC_GIGGLE_name <- NULL
for(cluster in as.numeric(colnames(BCC_GIGGLE))) 
{
	#' Use the transferred labels from scRNA-seq to annotate scATAC-seq clusters
	BCC_GIGGLE_name <- c(BCC_GIGGLE_name, names(sort(table(BCC_ATAC@meta.data[which(BCC_ATAC@meta.data$ATAC_snn_res.0.6==cluster),"assign.ident"]),decreasing=T)[1]))
}
colnames(BCC_GIGGLE) <- paste0(colnames(BCC_GIGGLE), ":", BCC_GIGGLE_name)
pdf("GSE129785_BCC_All_MAESTRO_ESRRA_TF_Enrichment_ATAC.pdf",width=8,height=5)
par(mar=c(10,5,3,3))
barplot(as.matrix(sort(BCC_GIGGLE['ESRRA',], decreasing = T))[1,], border=NA, col="black", las=2, ylab="GIGGLE enrichment score", main="Enrichment of ESRRA binding in BCC scATAC-seq clusters")
dev.off()

#' Scanning the enriched regulators in scRNA-seq clusters using LISA.
BCC_LISA <- read.delim('/homes/cwang/projects/DATA/SCRNAseq/Data_immunotherapy/GSE123814_human_aPD1/GSE123814_human_aPD1_bcc_TF_lisa.txt')
colnames(BCC_LISA) <- gsub("cluster","",colnames(BCC_LISA))
BCC_LISA_name <- NULL
for(cluster in as.numeric(colnames(BCC_LISA))) 
{
	#' Use the original labels from scRNA-seq to annotate the clusters
	BCC_LISA_name <- c(BCC_LISA_name, names(sort(table(BCC_RNA@meta.data[which(BCC_RNA@meta.data$RNA_snn_res.0.6==cluster),"assign.ident"]),decreasing=T)[1]))
}
colnames(BCC_LISA) <- paste0(colnames(BCC_LISA), ":", BCC_LISA_name)
pdf("GSE129785_BCC_All_MAESTRO_ESRRA_TF_Enrichment_RNA.pdf",width=8,height=5)
par(mar=c(10,5,3,3))
barplot(as.matrix(sort(BCC_LISA['ESRRA',], decreasing = T))[1,], border=NA, col="black", las=2, ylab="LISA enrichment score(-log10Pvalue)", main="Enrichment of ESRRA binding in BCC scRNA-seq clusters")
dev.off()















