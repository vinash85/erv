#NB data process
setwd('~/Dropbox/ML_input/NB_data/')
options(stringsAsFactors = F)
library(data.table)
library(tidyverse)
library(GEOquery)
library(sva)

#annotate the data 
NB_files = list.files(pattern = '.txt$')
for (nb.file in NB_files) {
  gse=getGEO(filename= nb.file)
  platform.info <- getGEO(gse@annotation, destdir=".")
  #the symbol name is not consistent between different platform 
  loc=grep(colnames(platform.info@dataTable@table),value = T,
           pattern = 'symbol',ignore.case = T)
  tmp=merge(platform.info@dataTable@table[,c('ID',loc)],
            gse@assayData$exprs,by.x=1,by.y=0)
  gene_list=tmp[,loc]
  tmp=tmp[,-c(1:2)]
  #get mean for each gene
  tmp.merged=aggregate(tmp,list(gene_id = gene_list),mean)
  filename=gsub(nb.file,pattern = '_series_matrix.txt',replacement = '')
  loc=which(tmp.merged$gene_id!='' & !is.na(tmp.merged$gene_id))
  tmp.merged=tmp.merged[loc,]
  write.table(tmp.merged,quote = F, row.names = F,sep = '\t',
              file =paste0(filename,'.annot.txt') )
}

#read + merge annot files and remove batch 
annot_files=list.files(pattern = '.annot')
#remove these platforms that the negative value cant be process by combat
#quantile normalized intensities GPL16876, which is a methods pf ranking columns(samples)
#others are samll sample size and high hetergenety
#part of GPL16876 are cell line, not human (649 neuroblastoma tumors)
human_sample=read.table('./GSE45547_acession.txt')
annot.data= fread(annot_files[2],data.table = F)
annot.data=annot.data[,colnames(annot.data) %in% c('gene_id',human_sample$V3)]
write.table(annot.data,quote = F, row.names = F,sep = '\t',
            file ='GSE45547_GPL16876_human.annot.txt' )

#TIDE input
#normalize the gene acrosss samples
rownames(annot.data) = annot.data$gene_id
annot.data = annot.data[,-1]
annot.data_nor = log2(annot.data)
annot.data_nor = annot.data_nor - rowMeans(annot.data_nor)
write.table(annot.data_nor,quote = F, row.names = T,sep = '\t',
            file ='GSE45547_GPL16876_TIDE_input.txt' )

#read cibersort and tide results
tide = read.csv('./TIDE_results_NB.csv')
tide_p = tide[,c(1,3,8:12)]
cibersort = fread('./CIBERSORT.Output_Job11.txt',data.table = F)
cibersort_tide= merge(cibersort,tide_p,by.x = 1,by.y=1)
write.table(cibersort_tide,quote = F, row.names = T,sep = '\t',
            file ='./Results/GSE45547_GPL16876_cibersort_tide.txt' )



#######
#!no use ; only choose the 649 tumors; dut to the combat error
for (i in 1:length(annot_files)) {
  annot.data= fread(annot_files[i],data.table = F)
  tmp.annot=data.frame(sample = colnames(annot.data)[-1],
                       batch=i)
  if (i==1) {
    Annot_Data = annot.data
    Sample.batch = tmp.annot
  }
  else{
    Annot_Data=merge(Annot_Data,annot.data,by='gene_id')
    Sample.batch=rbind(Sample.batch,tmp.annot)
  }
}
rownames(Annot_Data) = Annot_Data[,1]
Annot_Data = Annot_Data[,-1]
annot.data.comb = ComBat(data=as.matrix(Annot_Data),
                         batch = as.factor(Sample.batch$batch))
########

