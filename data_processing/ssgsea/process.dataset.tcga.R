
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
normalize.std = function(tt){
    (tt - min(tt, na.rm=T))/(max(tt,na.rm=T) - min(tt, na.rm=T))
}

process.dataset = function(dataset_ssgsea, pathway_order, dataset_phenotype, phenotype_order, output.dir, fix_patient_name =F, ICB_dataset =F) {

    library(data.table)
    dir.create(output.dir, showWarnings = FALSE)

    dataset_ssgsea = fread(dataset_ssgsea)
    pathway_order = fread(pathway_order)
    dataset_phenotype = fread(dataset_phenotype)
    xx = load(phenotype_order); phenotype_order = eval(parse(text=xx))
    # phenotype_order = fread(phenotype_order)
    dataset_ssgsea_mat= t(as.matrix(dataset_ssgsea[,seq(2,ncol(dataset_ssgsea)),with=F]))
    colnames(dataset_ssgsea_mat) = dataset_ssgsea$V1
# identical(dataset_ssgsea$V1, pathway_order$pathway)
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$pathway]  
    dataset_ssgsea_mat = dataset_ssgsea_mat[,pathway_order$order]
    patient.name = rownames(dataset_ssgsea_mat)
    patient.name = gsub(patient.name, pattern="-", replacement=".")
    # patient.name = substring(patient.name, 1, 12)
    rownames(dataset_ssgsea_mat) = patient.name

    ## there are duplicates in patient names because same patient have multiple expression. 

# phenotype data
    setnames(dataset_phenotype, 2, "patient.name")
    dataset_phenotype$patient.name = gsub(dataset_phenotype$patient.name, pattern="-", replacement=".")
    as.mynumeric = function(xx) as.numeric(ifelse(xx == '[Not Available]', NA, xx))
    # dataset_phenotype1 =dataset_phenotype
    cols = colnames(dataset_phenotype)[-(1:2)]
    dataset_phenotype[ , (cols) := lapply(.SD, as.mynumeric), .SDcols = cols ]
    only_in_phenotype = setdiff(dataset_phenotype$patient.name, patient.name)
    only_in_ssgsea = setdiff( patient.name, dataset_phenotype$patient.name)
    common.patients = intersect(patient.name, dataset_phenotype$patient.name)
    dataset_ssgsea_sel = dataset_ssgsea_mat[match(common.patients, patient.name), ] 

    phenotype_sel = dataset_phenotype[match(common.patients, dataset_phenotype$patient.name)]
# phenotype_sel[1:2,]

    colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern=" ", replacement="_")
    colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern="-", replacement="_")

    if(ICB_dataset){
        impute_reponse  = FALSE
        if(impute_reponse){
            ######
            # imputing  response
            #########
            genetech.patients = intersect(grep(phenotype_sel$patient.name, pattern="^SAM", value=T), rownames(dataset_ssgsea_sel))
            phenotype_sel.mod = phenotype_sel[match(genetech.patients, patient.name)]
            phenotype_sel.mod[, Response:=as.double(Response)]
            phenotype_sel.mod[is.na(Response) & (vital_status == 1) & (survive < 3)]$Response = 0
            phenotype_sel.mod[is.na(Response) & (survive > 7)]$Response = 1
            dataset_ssgsea_sel = dataset_ssgsea_sel[genetech.patients,]
            phenotype_mat =  as.matrix(phenotype_sel.mod[,2:4,with=F])
            phenotype.ext.mat = phenotype_mat 
            

            phenotype_sel.mat = as.matrix(phenotype_sel.mod)
            phenotype_mat = phenotype_sel.mat[,match(phenotype_order[seq(length(phenotype_order) -1)], colnames(phenotype_sel.mat)) ]
            phenotype.ext.mat = apply(as.matrix(cbind(phenotype_mat, phenotype_sel.mod$Response)),2,as.numeric)




            }else{

                # create groups
                setnames(phenotype_sel, "survive", "survive_ICB")
                setnames(phenotype_sel, "vital_status", "vital_status_ICB")
                prefix = c("^H_VB", "^SRR", "^ERR", "^CA", "^PD", "^SAM")
                phenotype_sel$dataset= rep(NA, nrow(phenotype_sel))
                for (pre in seq(length(prefix))) {
                    pre.curr = prefix[pre]
                    inx.curr = grep(pre.curr, phenotype_sel$patient.name)
                    phenotype_sel$survive_ICB[inx.curr] = normalize.std(phenotype_sel$survive_ICB[inx.curr])
                    phenotype_sel$dataset[inx.curr] = pre


                }
# O is non-responder and 1 is responder 
                phenotype_sel[,response:= ifelse(Response %in% c("1", "CR", "PR"), 1, ifelse(Response %in% c("0", "SD", "PD"),0, NA))]

                response.df = phenotype_sel[,list(survive_ICB, vital_status_ICB,  response)] 

                phenotype_mat =  as.matrix(phenotype_sel[,-1,with=F])
                phenotype_mat = phenotype_mat[,match(phenotype_order[seq(length(phenotype_order) -3)], colnames(phenotype_mat)) ]
                phenotype.ext.mat = apply(as.matrix(cbind(phenotype_mat, response.df)),2,as.numeric)

            }
        }
        



        }else{

    # setnames(phenotype_sel, "OS", "survive")
    # setnames(phenotype_sel, "OS.Event", "vital_status")

            normalize.survival = T
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
        }


        write.table(file=paste0(output.dir, "/dataset_ssgsea.txt"),x = dataset_ssgsea_sel,
            row.names = F, col.names =T,  sep="\t", quote=F )
        write.table(file=paste0(output.dir, "/dataset_phenotype.txt"),x = phenotype.ext.mat,
            row.names = F, col.names =T,  sep="\t", quote=F )



        rand_inx = sample(nrow(dataset_ssgsea_sel))
        dataset_ssgsea_sel_shuffle = dataset_ssgsea_sel[rand_inx,]
        phenotype.ext.mat_shuffle = phenotype.ext.mat[rand_inx,]
        train.inx = 1:ceiling(.7 * nrow(dataset_ssgsea_sel_shuffle))
        val.inx = ceiling(.7 * nrow(dataset_ssgsea_sel_shuffle)): ceiling(.9 * nrow(dataset_ssgsea_sel_shuffle))
        test.inx = ceiling(.9 * nrow(dataset_ssgsea_sel_shuffle)):nrow(dataset_ssgsea_sel_shuffle)
        
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

        write.table(file=paste0(output.dir, "/sgsea_test.txt"),x = dataset_ssgsea_sel_shuffle[test.inx,],
            row.names = F, col.names =T,  sep="\t", quote=F )
        write.table(file=paste0(output.dir, "/phenotype_test.txt"),x = phenotype.ext.mat_shuffle[test.inx,],
            row.names = F, col.names =T,  sep="\t", quote=F )

    }

if(FALSE){
    phenotype_new_order = c("cancertype", phenotype_order[1:33], "oxphos_score",phenotype_order[34:41] )
    phenotype_order = phenotype_new_order
    save(file="/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData", phenotype_order)

}