
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
rownames(dataset_ssgsea_mat) = patient.name

# phenotype data
setnames(dataset_phenotype, 1, "patient.name")
dataset_phenotype$patient.name = gsub(dataset_phenotype$patient.name, pattern="-", replacement=".")
only_in_phenotype = setdiff(dataset_phenotype$patient.name, patient.name)
only_in_ssgsea = setdiff( patient.name, dataset_phenotype$patient.name)
common.patients = intersect(patient.name, dataset_phenotype$patient.name)
dataset_ssgsea_sel = dataset_ssgsea_mat[match(common.patients, patient.name), ] 

phenotype_sel = dataset_phenotype[match(common.patients, dataset_phenotype$patient.name)]
# phenotype_sel[1:2,]

colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern=" ", replacement="_")
colnames(phenotype_sel) = gsub(colnames(phenotype_sel), pattern="-", replacement="_")

if(ICB_dataset){
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



}else{

    # setnames(phenotype_sel, "OS", "survive")
    # setnames(phenotype_sel, "OS.Event", "vital_status")
    phenotype_mat =  as.matrix(phenotype_sel[,-1,with=F])

    phenotype.ext.mat = phenotype_mat[,match(phenotype_order, colnames(phenotype_mat)) ]
}

                        


}

