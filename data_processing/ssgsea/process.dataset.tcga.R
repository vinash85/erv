
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
source("./source.R")

process.dataset = function(dataset_ssgsea, pathway_order, dataset_phenotype, phenotype_order, output.dir, fix_patient_name =F, ICB_dataset =F) {

    library(data.table)
    dir.create(output.dir, showWarnings = FALSE)

    dataset_ssgsea = fread(dataset_ssgsea)
    pathway_order = fread(pathway_order)
    dataset_phenotype = fread(dataset_phenotype)
    xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
    # phenotype_order = fread(phenotype_order)
    dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
    # setnames(dataset_ssgsea, 1, "gene_name")
# identical(dataset_ssgsea$V1, pathway_order$pathway)
  
    tpm = T
    if(!tpm){
    colnames(dataset_ssgsea_mat) = dataset_ssgsea$V1
        dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$pathway]  
        dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$order]
        }else{
            colnames(dataset_ssgsea_mat) = dataset_ssgsea$gene_name
            pcgs = fread("/liulab/asahu/data/ssgsea/xiaoman/./pcg.txt")
           # dataset_ssgsea_mat = dataset_ssgsea_mat[,toupper(colnames(dataset_ssgsea_mat)) %in% toupper(pcgs$Gene)] 
            load("/liulab/asahu/data/ssgsea/xiaoman/commmon.genes.RData")
            dataset_ssgsea_mat = dataset_ssgsea_mat[ ,common.genes] 
            stopifnot(any(!is.na(dataset_ssgsea_mat)))

        }


 
    patient.name = rownames(dataset_ssgsea_mat)
    patient.name = gsub(patient.name, pattern="-", replacement=".")
    rownames(dataset_ssgsea_mat) = patient.name

    setnames(dataset_phenotype, 2, "patient.name")
    dataset_phenotype$patient.name = gsub(dataset_phenotype$patient.name, pattern="-", replacement=".")
    as.mynumeric = function(xx) as.numeric(ifelse(xx == '[Not Available]', NA, xx))
    cols = setdiff(colnames(dataset_phenotype)[-(1:2)], "cancertype")
    dataset_phenotype[ , (cols) := lapply(.SD, as.mynumeric), .SDcols = cols ]
    only_in_phenotype = setdiff(dataset_phenotype$patient.name, patient.name)
    only_in_ssgsea = setdiff( patient.name, dataset_phenotype$patient.name)
    common.patients = intersect(patient.name, dataset_phenotype$patient.name)
    dataset_ssgsea_sel = dataset_ssgsea_mat[match(common.patients, patient.name), ] 

    phenotype_sel = dataset_phenotype[match(common.patients, dataset_phenotype$patient.name)]

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
            # phenotype.ext.mat = phenotype_mat[,match(phenotype_order, colnames(phenotype_mat)) ]
            phenotype.ext.mat = phenotype.ext.mat[,match(phenotype_order, colnames(phenotype.ext.mat)),with=F ]

        if(pca){
            temp_out = get_pca(dataset_ssgsea_sel, subsample=.2) 
            # std_cs = cumsum(temp_out$pca_obj$sdev^2/sum(temp_out$pca_obj$sdev^2))
            # sum(std_cs < 0.95)
            # sum(std_cs < 0.9)
            # sum(std_cs < 0.93)
            pca_obj = temp_out$pca_obj
            pca_obj$len_selected = 50
            save(file=paste0(output.dir, "/pca_obj.RData"), pca_obj)
            pca_out_sel = temp_out$pca_out[,seq(pca_obj$len_selected)]
            dataset_ssgsea_sel = pca_out_sel 
        }
        if(combine_phenotype_feature){
            barcode = substring(rownames(dataset_ssgsea_sel), 1,12)
            # library(readxl)
            # phenotype = excel_read("liulab/asahu/data/ssgsea/xiaoman/mmc2.xlsx")
            phenotype = fread("/liulab/asahu/data/ssgsea/xiaoman/mmc2.txt")
            cols = colnames(phenotype)
            cols[duplicated(cols)] = paste0(cols[duplicated(cols)], ".1")
            cols = gsub(cols, pattern=" ", replacement=".")
            colnames(phenotype) = cols

            phenotype$TCGA.Participant.Barcode = gsub(phenotype$TCGA.Participant.Barcode, pattern="-", replacement=".")
            sum(barcode %in% phenotype$TCGA.Participant.Barcode )
            reorder = match(barcode, phenotype$TCGA.Participant.Barcode)
            feature.sel = c(3:36)
            f1 = phenotype[reorder, feature.sel, with=F]
            new.pheno = cbind(phenotype.ext.mat, f1)
            phenotype.ext.mat = new.pheno

        }

        write.table(file=paste0(output.dir, "/dataset_ssgsea.txt"),x = dataset_ssgsea_sel,
            row.names = F, col.names =T,  sep="\t", quote=F )
        write.table(file=paste0(output.dir, "/dataset_phenotype.txt"),x = phenotype.ext.mat,
            row.names = F, col.names =T,  sep="\t", quote=F )



        rand_inx = sample(nrow(dataset_ssgsea_sel))
        dataset_ssgsea_sel_shuffle = dataset_ssgsea_sel[rand_inx,]
        phenotype.ext.mat_shuffle = phenotype.ext.mat[rand_inx,]
        train.inx = 1:ceiling(.85 * nrow(dataset_ssgsea_sel_shuffle))
        # val.inx = ceiling(.8 * nrow(dataset_ssgsea_sel_shuffle)): ceiling(.9 * nrow(dataset_ssgsea_sel_shuffle))
        val.inx = ceiling(.85 * nrow(dataset_ssgsea_sel_shuffle)):nrow(dataset_ssgsea_sel_shuffle)
        
        # train.inx = 1:ceiling(.5 * nrow(dataset_ssgsea_sel))
        # val.inx = ceiling(.5 * nrow(dataset_ssgsea_sel)): ceiling(nrow(dataset_ssgsea_sel))

        write.table(file=paste0(output.dir, "/ssgsea_train.txt"),x = dataset_ssgsea_sel_shuffle[train.inx,],
            row.names = F, col.names =T,  sep="\t", quote=F )
        write.table(file=paste0(output.dir, "/phenotype_train.txt"),x = phenotype.ext.mat_shuffle[train.inx,],
            row.names = F, col.names =T,  sep="\t", quote=F )

        write.table(file=paste0(output.dir, "/ssgsea_val.txt"),x = dataset_ssgsea_sel_shuffle[val.inx,],
            row.names = F, col.names =T,  sep="\t", quote=F )
        write.table(file=paste0(output.dir, "/phenotype_val.txt"),x = phenotype.ext.mat_shuffle[val.inx,],
            row.names = F, col.names =T,  sep="\t", quote=F )

        # write.table(file=paste0(output.dir, "/sgsea_test.txt"),x = dataset_ssgsea_sel_shuffle[test.inx,],
        #     row.names = F, col.names =T,  sep="\t", quote=F )
        # write.table(file=paste0(output.dir, "/phenotype_test.txt"),x = phenotype.ext.mat_shuffle[test.inx,],row.names = F, col.names =T,  sep="\t", quote=F )

    }

if(FALSE){
    phenotype_new_order = c("cancertype", phenotype_order[1:33], "oxphos_score",phenotype_order[34:41] )
    phenotype_order = phenotype_new_order
    save(file="/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData", phenotype_order)

}