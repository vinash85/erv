# CUDA_VISIBLE_DEVICES=3 python train.py  --data_dir  ~/project/deeplearning/icb/data/tcga/tf/binarize/datasets_tsne_list.txt --model_dir ~/project/deeplearning/icb/data/tcga/tf/binarize/.  --tensorboard_prefix surv_noaug_bin_ --hyper_param "params.add_noise=0"
# 
# 
col.names = colnames(fread("~/project/deeplearning/icb/data/tcga/tf/binarize/dataset_val.txt",nrows = 1))
start = 4001 + 244 -2
pheno.all.inx = 4001 + 244:438
pheno.type = c(rep("continuous", 400-244 +1), rep("binary",430-401+1), rep( "survival",438-431+1))
pheno.next = start+1
pheno.last.index = pheno.all.inx[length(pheno.all.inx)]
commands = list()
cuda.id = 0
while(pheno.next < pheno.last.index){
    label.curr = paste0("Pheno_",pheno.next, "_", col.names[pheno.next], "_")
    if (pheno.type[pheno.next - start] == "continuous"){
        pheno.str=sprintf("continuous_phenotype_indices=[%d]", pheno.next-2)
        pheno.next = pheno.next + 1
    }else if (pheno.type[pheno.next - start] == "binary"){
        pheno.str=sprintf("binary_phenotype_indices=[%d]", pheno.next-2)
        pheno.next = pheno.next + 1
    }else if(pheno.type[pheno.next - start]=="survival"){
        pheno.str=sprintf("survival_indices=[%d, %d]", pheno.next-2, pheno.next-2+1)
        pheno.next = pheno.next + 2
    }
    commands[[label.curr]] = sprintf("CUDA_VISIBLE_DEVICES=%d python train.py  --data_dir  ~/project/deeplearning/icb/data/tcga/tf/binarize/datasets_tsne_list.txt --model_dir ~/project/deeplearning/icb/data/tcga/tf/binarize/.  --tensorboard_prefix %s --hyper_param \"params.%s\"",cuda.id, label.curr, pheno.str)
    cuda.id = (cuda.id +1) %% 4
    
}


library(parallel)
# Using fork()
mc.cores = 16 
out = list()
for (ii in seq(ceiling(length(commands)/mc.cores))) {
    last.task.in.itr = min(length(commands), mc.cores * (ii))
    tasks.curr = commands[(mc.cores * (ii-1) +1 ): last.task.in.itr] 
    print(sprintf("running commands %d", last.task.in.itr))
    # print(tasks.curr)
out[[ii]] <- mclapply(tasks.curr, function(tt) system(c(unlist(tt)),intern =  T),
    mc.cores = mc.cores)
# out <- mclapply(tasks.curr, function(tt) c(unlist(tt)),
                # mc.cores = mc.cores, mc.cleanup=TRUE)
}
