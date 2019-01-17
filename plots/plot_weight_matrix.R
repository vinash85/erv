library(data.table)
library(ggplot2)

model_dir = "./experiments/reproduce_tcga_43_l1/"
plot.dir = paste0(model_dir, "plots")
load("/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData")

dir.create(plot.dir)
weight.file = paste0(model_dir, "linear_output_weights.txt")
weights = as.matrix(fread(weight.file))
rownames(weights) = c( "survival", tcga.phenotypes[seq(3,nrow(weights)+1)])
colnames(weights) = paste0("e", seq(1,ncol(weights)))

dat = t(weights)
# colnames(dat) = tcga.phenotypes
sds = apply(dat,1, sd)
dat = dat[sds > 0, ]
hc = hclust(as.dist(1-cor(dat, method="spearman")), method="complete")
hr = hclust(as.dist(1-cor(t(dat), method="spearman")), method="complete")


library(RColorBrewer)
# pats.new.col = pats.new; levels(pats.new.col) = brewer.pal(length(levels(pats.new)),'Set1')
# site.new.col = site.new; levels(site.new.col) = brewer.pal(length(levels(site.new)),'Set2')
# type.new.col = type.new; levels(type.new.col) = brewer.pal(length(levels(type.new)),'Set3')
# ColSideColors<-cbind( patients=as.character(pats.new.col), type=as.character(type.new.col), site=as.character(site.new.col))
library(heatmap3)

pdf(paste0(plot.dir,"/output_linear_weights.pdf"), width = 12, height = 12)
heatmap3(dat, Rowv=as.dendrogram(hr),  Colv=as.dendrogram(hc), scale="row", balanceColor=T, showRowDendro=F, labRow=NA, ColSideCut=0.6, margins=c(8,8) )
dev.off()


library(corrplot)


weights_non_empty = weights[!grepl(pattern="empty", rownames(weights)),]


dat = weights_non_empty
# dat = dat[,order(new.labs)]
hc = hclust(as.dist(1-cor(dat, method="spearman")), method="complete")
dat = dat[,hc$order]
M<-cor(dat,method="spearman")

# col <- c( rep("#FFFFFF", 25), rep("gray67",9), rep("khaki", 1),  "yellow", "yellow", "orange", "red", "red4")
# pdf("ssgsea_correlation_combat_categorical.pdf")
pdf(paste0(plot.dir,"/embedding_correlation.pdf"), width = 12, height = 12)
corrplot(M, 
	# col=col,
         type="upper",
          # order="hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         tl.cex = 0.5,
         number.cex = .5,
         # Combine with significance
         # p.mat = p.mat,
         # sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
         )
dev.off()



dat = t(weights_non_empty)
hc = hclust(as.dist(1-cor(dat, method="spearman")), method="complete")
dat = dat[,hc$order]
M<-cor(dat,method="spearman")

# col <- c( rep("#FFFFFF", 25), rep("gray67",9), rep("khaki", 1),  "yellow", "yellow", "orange", "red", "red4")
# pdf("ssgsea_correlation_combat_categorical.pdf")
pdf(paste0(plot.dir,"/output_correlation.pdf"), width = 12, height = 12)
corrplot(M, 
	# col=col,
         type="upper",
          # order="hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         tl.cex = 0.5,
         number.cex = .5,
         # Combine with significance
         # p.mat = p.mat,
         # sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
         )
dev.off()

