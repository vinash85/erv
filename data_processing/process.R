setwd('d:/Harvard/Lab/Avi/OXPHOS/')
options(stringsAsFactors = F)
library(data.table)
clinicalData=read.table('d:/Harvard/Lab/Avi/OXPHOS/IMvigor210.clinical',row.names = 1,header = T,sep = '\t')
colnames(clinicalData)[1:2]=c('survive','vital_status')
#from RSEM TPM
GeneExpr=read.delim('d:/Harvard/Lab/Avi/OXPHOS/Genetech_expression_TPM.txt')

####################################
#cluster ssgsea

myC7<-getGmt("d:/Harvard/Lab/Avi/RNAseq/Braf/c7.all.v6.2.symbols.gmt", geneIdType=SymbolIdentifier(),collectionType=BroadCollection(category="c7"), sep="\t")

esrnaseq.genetech <- gsva(as.matrix(GeneExpr), myC7 , min.sz=5, max.sz=500,
                          kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=2, ssgsea.norm=TRUE)
esrnaseq.genetech=as.data.frame(esrnaseq.genetech)
#need scale?
#using Absolute distance to cluster pathways
dist.genetech=dist(esrnaseq.genetech,method = 'manhattan')
hclust_genetech <- hclust(dist.genetech, method = 'complete')
gsea.out=data.frame(pathway=hclust_genetech$labels, order=hclust_genetech$order)
write.table(gsea.out,file = './others/ssgsea.order_Genetech.txt',row.names = F,quote = F,sep = '\t')
write.table(esrnaseq.genetech,file = './Avi/ssgsea_Genetech.txt',row.names = F,quote = F,sep = '\t')

######################################
#biomarkers:

pd1=as.data.frame(t(GeneExpr[rownames(GeneExpr)%in%c("PDCD1","CD274"),]))
ciber.genetech=read.csv('d:/Harvard/Lab/Avi/OXPHOS/CIBERSORT.Output_Genetech.csv',row.names = 1)
ciber.genetech=ciber.genetech[rownames(pd1),]
tide.genetech=read.csv('d:/Harvard/Lab/Avi/OXPHOS/cohort/ICB_for_Jingxin/TIDE/TIDE_genetech.csv',row.names = 1)
tide.genetech=tide.genetech[rownames(pd1),]

biom.genetech=clinicalData[,c('survive','vital_status','Response','Mutation')]
biom.genetech=cbind(biom.genetech,pd1,MSI=rep(NA,nrow(pd1)),t(esrnaseq.genetech),ciber.genetech,tide.genetech)
write.table(biom.genetech,file = './others/biomarker_Genetech.txt',quote = F,sep = '\t')

############################################################################
#for tcga
#note merge by biopsy id not the patient id
#firstly merge respective id, and then last step to merge biopsy id and patient id 
tcga.clinical=fread('d:/Harvard/Lab/Peter/clinical_PANCAN_patient_with_followup.tsv',data.table = F)

load('d:/Harvard/Lab/Avi/RNAseq/Braf/data/clinicaldata_sub.RData')
biom.tcga=cbind(Clinicaldata_sub,Response=rep(NA,nrow(Clinicaldata_sub)),MSI=tcga.clinical$microsatellite_instability)

setwd('d:/Harvard/Lab/Avi/OXPHOS/others/Results/Tumor_Dysf_Excl_scores/')
Timerfiles=list.files('d:/Harvard/Lab/Avi/OXPHOS/others/Results/Tumor_Dysf_Excl_scores/',
                      pattern = 'RNASeq.norm_subtract.OS_base')
timerData=lapply(Timerfiles, function(x){data.frame(fread(x,data.table = F), row.names=1)})
timer=do.call(rbind,timerData)
rownames(timer)=gsub(rownames(timer),pattern = '-',replacement = '.')
tcga.biom=merge(biom.tcga,timer,by.x=1,by.y=0,all=T)

#biospy id
#biomarkers:
#mutation 0
tcga.mutation=fread('d:/Harvard/Lab/Avi/OXPHOS/others/mc3.v0.2.8.PUBLIC.nonsilentGene.xena',data.table=F)
#colnames(tcga.mutation)=substr(colnames(tcga.mutation),1,12)
ml.tcga=data.frame(Mutation=colSums(tcga.mutation[,-1]))
rownames(ml.tcga)=colnames(tcga.mutation[,-1])
rownames(ml.tcga)=gsub(rownames(ml.tcga),pattern = '-',replacement = '.')
#cibersort 0
setwd('d:/Harvard/Lab/Avi/OXPHOS/others/TCGA_cibersort_out/')
ciberfiles=list.files(path = 'd:/Harvard/Lab/Avi/OXPHOS/others/TCGA_cibersort_out/')
ciberData=lapply(ciberfiles, function(x){data.frame(fread(x,data.table = F), row.names=1)})
ciber=do.call(cbind,ciberData)
ciber=as.data.frame(t(ciber))
#pd1 0
pd1=as.data.frame(t(read.table('../pd_expr.txt',header = T,row.names = 1)))
#cancertype
#oxphos level
oxphos_level=read.csv('d:/Harvard/Lab/Avi/OXPHOS/loading_data/TIDE_Oxphos/oxphosLevel_allcancers.csv')
oxphos_level$cancer=gsub(oxphos_level$cancer,pattern = 'TCGA_',replacement = '')
colnames(oxphos_level)[c(3,7)]=c('oxphos_score','cancertype')
rownames(oxphos_level)=oxphos_level$object
oxphos_level=oxphos_level[,c('oxphos_score','cancertype')]
#merge all biopsy id
#tcga.biom.biopsy=Reduce(function(x,y) merge(x=x,y=y,by=0),
#                        list(ciber,pd1,oxphos_level))
tmp=merge(ciber,pd1,by=0)
tmp2=merge(tmp,oxphos_level,by.y=0,by.x = 1)
tmp2$Row.names=substr(tmp2$Row.names,1,15)
tcga.biom.biopsy=merge(tmp2,ml.tcga,by.x=1,by.y=0,all.x = T)
tcga.biom.biopsy$patientID=substr(tcga.biom.biopsy$Row.names,1,12)
tcga.biom.final=merge(tcga.biom.biopsy,tcga.biom,by.y = 1,by.x='patientID',all.x = T)
colnames(tcga.biom.final)[2]='BiopsyID'
tcga.biom.final=tcga.biom.final[,c(1:2,33:36,32,3:29,37:41,31,30)]
write.table(tcga.biom.final,file = 'd:/Harvard/Lab/Avi/OXPHOS/others/tcga_biom_oxphos.txt',
            row.names = F,quote = F,sep = '\t')
save.image('tcga_biom_ssgsea.RData')

################
#ssgsea
#cluster all cancers
setwd('d:/Harvard/Lab/Avi/OXPHOS/others/')
#setwd('/home/xw165/ssgsea/')
gseafile=list.files(pattern = 'ssgseaMat')

options(stringsAsFactors = F)
gesaData=lapply(gseafile, function(x){data.frame(fread(x,header = F,data.table = F), row.names=1)})
gsea=do.call(cbind,gesaData)

#using Absolute distance to cluster pathways
dist.panc=dist(gsea,method = 'manhattan')
hclust_panc <- hclust(dist.panc, method = 'complete')
gsea.out=data.frame(pathway=hclust_panc$labels, order=hclust_panc$order)
write.table(gsea.out,file = 'ssgsea.order_panc.txt',row.names = F,quote = F,sep = '\t')
  