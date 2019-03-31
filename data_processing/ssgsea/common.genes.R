#### TCGA BLCA + Neo-antigen prediction. 
# find common genes in genentech and tcga 
aa = fread("/liulab/asahu/data/ssgsea/xiaoman/Genetech_genes.txt")
bb = fread("/liulab/asahu/data/ssgsea/xiaoman/TCGA_genes.txt")
cc = fread("/liulab/asahu/data/ssgsea/xiaoman/pcg.txt")

common.genes = intersect(intersect(unlist(aa), unlist(bb)), cc$Gene)

save(common.genes, file="/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")