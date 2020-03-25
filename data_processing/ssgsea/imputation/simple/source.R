#' Get list of selected immune genes from literature
#'
#' @return
#' @export
#'
#' @examples
get.immune.genes = function(){
    T_Cell_extra = c("IL10", "IDO", "TGFB1", "TGFB2", "TGFBR1", "TGFBR1", "CD37", "TLR", "Arginase")
    APC_2 = c("A2AR", "VISTA", "B7_h3", "PDL1", "PDL2", "CD80", "CD86", "Galectin_9", "Ox40L", "CD40", "B7RP1", "CD70", "HVEM", "GITRL", "TNFSF9", "CD155", "CD112")
    T_Cell_1 = c("CTLA4", "TIM3", "OX40", "CD40L", "ICOS", "CD27", "BTLA", "LAG3", "TCR", "KIR", "GITR", "TNFRSF9", "CD226", "TIGIT")
    checkpoint.genes = unique(c(T_Cell_extra, APC_2, T_Cell_1))
    load("/liulab/asahu/data/ssgsea/xiaoman/getz/all.tcga.genes.RData")
    setdiff(checkpoint.genes, all.genes)
    checkpoint.genes.1 = intersect(checkpoint.genes, all.genes)
    
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
    
    checkpoint.genes.semifinal = unique(c(checkpoint.genes.1, checkpoint.genes.rescue, pd1.genes, Galectin_9.genes, stat.genes, TNF.genes, il2.gene, il7.gene, il4.gene, il6.gene, il10.gene, HAVC.gene, gzm.genes, traf.genes, nfk.genes, cd40.genes, igh.genes, cd3.genes,tra.genes, kir.genes, tgf.genes, antigen.presentation.genes, serpin.genes, vegf.genes))
    immune.genes = unique(intersect(checkpoint.genes.semifinal, all.genes))    
    immune.genes
}

#' Get immune landscape and msi data
#'
#' @param barcode = barcode of TCGA samples 
#'
#' @return
#' @export
#'
#' @examples
get.tcga.immune.landscape <- function(barcode) {
    ## neoantigen 
    phenotype = fread("/liulab/asahu/data/ssgsea/xiaoman/mmc2.txt")
    cols = colnames(phenotype)
    cols[duplicated(cols)] = paste0(cols[duplicated(cols)], ".1")
    cols = gsub(cols, pattern=" ", replacement=".")
    colnames(phenotype) = cols
    
    phenotype$TCGA.Participant.Barcode = gsub(phenotype$TCGA.Participant.Barcode, pattern="-", replacement=".")
    # sum(barcode %in% phenotype$TCGA.Participant.Barcode )
    reorder = match(barcode, phenotype$TCGA.Participant.Barcode)
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
    
    # length(sort(table(phenotype.ext.mat$cancertype)))
    # cancertype.dt = data.table(ID = rownames(dataset_ssgsea_sel), cancertype = as.factor(phenotype.ext.mat$cancertype))
    # cancertype.pheno = mltools::one_hot(cancertype.dt, sparsifyNAs =T)
    # 
    # cancertype.pheno = cancertype.pheno[,-1, with=F]
    
    ## select 10 genes per phenotype and take PC
    msi.neoantigen = cbind(neoantigen.pheno, msi.pheno[,c("Total_nb_MSI_events", "MSI_exonic"), with=F]) 
    
    cors = cor(msi.neoantigen, dataset_ssgsea_sel.back,  use = "pairwise.complete.obs")
    top.cors = lapply(seq(nrow(cors)), function(tt) {
        xx = cors[tt,] 
        xx[order( abs(xx), decreasing = T)[1:20]] 
    }
    )
    names(top.cors)  = rownames(cors)
    # save(file="/liulab/asahu/data/ssgsea/xiaoman/ICB_datasets/top_correlated_genes_with_ICB_biomarkers.RData", top.cors)
    
    # save(file="/liulab/asahu/data/ssgsea/xiaoman/ICB_datasets/correlated_genes_with_ICB_biomarkers.RData", cors)
    
    top.genes = unique(unlist(lapply(top.cors, names)))
    list(neoantigen.pheno = neoantigen.pheno, msi.pheno = msi.pheno, top.correlated.genes = top.genes)
}

create.rrpa.dataset = function(expression_match, nontranscriptome.immune.factors, tcga.immune.factors){
    immune.genes = get.immune.genes()
    # nontranscriptome.immune.factors = tcga.dataset[,c(102:103, 110:123, 128:137),with=F]                                     
    # tcga.immune.factors = tcga.dataset[,c(66:137, 178:188),with=F]
    tcga.immune.genes = expression_match[,immune.genes]
    
    rrpa = fread("/liulab/asahu/data/ssgsea/xiaoman/TCGA_gdc/TCGA-RPPA-pancan-clean.txt")
    rrpa[,SampleID:=gsub(substring(SampleID,1,16), pattern="-", replacement = ".")]
    rrpa.match = rrpa[match(sample.name$x, SampleID)]
    rrpa.mat = rrpa.match[,-(1:2),with=F]
    temp = calc.cor(rrpa.mat, nontranscriptome.immune.factors, use = "pairwise.complete.obs")
    temp2 = calc.cor(rrpa.mat, tcga.immune.factors, use = "pairwise.complete.obs")
    temp3 = calc.cor(rrpa.mat, tcga.immune.genes, use = "pairwise.complete.obs")
    protien.sel1 = which(rowSums(abs(temp$estimate)> 0.3 , na.rm = T) > 0)
    protien.sel2 = which(rowSums(abs(temp2$estimate)> 0.4 , na.rm = T) > 0)
    protien.sel3 = which(rowSums(abs(temp2$estimate)> 0.4 , na.rm = T) > 0)
    protien.sel = unique(c(names(protien.sel1), names(protien.sel2), names(protien.sel3)))
    
    
    rrpa.dataset = rrpa.mat[,protien.sel,with=F]
    setnames(rrpa.dataset, seq(ncol(rrpa.dataset)), paste(colnames(rrpa.dataset), "prot", sep=".")) 
    
    cors.list = cor(rrpa.dataset, expression_match, use = "pairwise.complete.obs")
    cors.list = cors.list[colnames(rrpa.dataset),]
    aaa = data.table(rownames(cors.list), seq(nrow(cors.list)))
    topmost.cors.genes = sapply(aaa$V1, function(tt) {
        xx = cors.list[tt,] 
        xx = xx[abs(xx) > 0.15]
        names(which.max(abs(xx)))
    }
    )
    
    top.cors = lapply(aaa$V1, function(tt) {
        xx = cors.list[tt,] 
        xx = xx[order( abs(xx), decreasing = T)[1:10]] 
        xx[abs(xx) >0.15]
        names(xx)
    }
    )
    rrpa.topmost.cors = unique(unlist(topmost.cors.genes))
    rrpa.top.cors = unique(unlist(top.cors))
    out = list(rrpa.dataset=rrpa.dataset, rrpa.topmost.cors=rrpa.topmost.cors, rrpa.top.cors=rrpa.top.cors)
    out
}

get.tcga.erv = function(common.patients){
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
}


get.tcga.survival <- function(common.patients) {
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
    clinical.final.dt
}


get.tcga.mutation <- function(common.patients) {
    sample.name = list(x=common.patients)
    mut.genes = c("STK11", "B2M", "PTEN", "PTPN2", "APLNR")
    setdiff(mut.genes, all.genes)
    baf.complex.genes=fread("/liulab/asahu/data/ssgsea/xiaoman/baf_complex.csv")$V2
    baf.complex.genes =baf.complex.genes[-1]
    all.mut.genes = unique(c(mut.genes, baf.complex.genes))
    
    mutation = fread("/liulab/asahu/data/ssgsea/xiaoman/TCGA_gdc/TCGA_consolidated.abs_mafs_truncated.fixed.txt")
    head(mutation$Hugo_Symbol)
    # names(table(mutation.dt$Consequence))
    synonymous.mut.type = c("downstream_gene_variant",   "intron_variant",       "stop_retained_variant",    "synonymous_variant",       "upstream_gene_variant")   
    
    mutation[,SampleID:=gsub(substring(sample,1,16), pattern="-", replacement = ".")]
    mut.samples = unique(mutation$SampleID)
    mutation.dt = mutation[Hugo_Symbol %in% all.mut.genes]
    mutation.dt = mutation.dt[!(Consequence %in% synonymous.mut.type)]
    mutation.dt =mutation.dt[SampleID%in%sample.name$x]
    mut.bin.mat = table(mutation.dt[,.( SampleID , Hugo_Symbol)])
    mut.ccf.mat = xtabs(ccf_CI95_low~ SampleID + Hugo_Symbol, data=mutation.dt)
    mut.cff.mat = mut.ccf.mat[,colnames(mut.bin.mat)]
    
    mut.bin = matrix(NA, nrow=length(sample.name$x), ncol = ncol(mut.bin.mat), dimnames = list(sample.name$x, colnames(mut.bin.mat)))
    mut.bin[intersect(mut.samples, sample.name$x),] =0
    mut.bin[rownames(mut.bin.mat),] = mut.bin.mat
    
    mut.ccf = matrix(NA, nrow=length(sample.name$x), ncol = ncol(mut.ccf.mat), dimnames = list(sample.name$x, colnames(mut.ccf.mat)))
    mut.ccf[intersect(mut.samples, sample.name$x),] =0
    mut.ccf[rownames(mut.ccf.mat),] = mut.ccf.mat
    colnames(mut.ccf) = paste0(colnames(mut.ccf.mat), ".ccf")
    mut.dataset = cbind( mut.ccf, mut.bin)
    genes.curr1 = names(sort(colSums(mut.bin.mat),decreasing = T))[1:6]
    genes.curr = unique(c(mut.genes, genes.curr1))
    
    cors.list = cor(mut.dataset, expression_match, use = "pairwise.complete.obs")
    rownames(cors.list) = colnames(mut.dataset)
    aaa = data.table(rownames(cors.list), seq(nrow(cors.list)))
    aaa = aaa[V1 %in% c(genes.curr, paste0(genes.curr, ".ccf"))]
    
    
    topmost.cors.genes = sapply(aaa$V1, function(tt) {
        xx = cors.list[tt,] 
        xx = xx[abs(xx) > 0.1]
        names(which.max(abs(xx)))
    }
    )
    
    top.cors = lapply(aaa$V1, function(tt) {
        xx = cors.list[tt,] 
        xx = xx[order( abs(xx), decreasing = T)[1:10]] 
        xx[abs(xx) >0.1]
        names(xx)
    }
    )
    mut.topmost.cors = unique(unlist(topmost.cors.genes))
    mut.top.cors = unique(unlist(top.cors))
    list(mut.dataset, mut.top.cors, mut.topmost.cors)
}
