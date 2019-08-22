###
#######
#processing of datasets
######

## TCGA dataset 
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/TCGA_ssgsva.txt"
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/TCGA_ALLTPM.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/tcga_biom_oxphos.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/tcga/PCA/"

output.dir = "~/project/deeplearning/icb/data/tcga.blca/neoantigen/"

# ICB dataset

dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Avin/ICB_GSVA.txt"
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Genetech_expression_TPM.txt"


pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_ICB_oxphos.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/all_icb"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/pca/"
output.dir = "~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/"

dataset.prefix = "Mel_CTLA4_Chiappinelli"
dataset.prefix = "mGC_PD1_Kim"
dataset.prefix = "RCC_PD1_Miao"
# dataset_ssgsea = sprintf("/liulab/asahu/data/ssgsea/xiaoman/expression/%s.expression",dataset.prefix)
dataset_ssgsea = sprintf("/liulab/asahu/data/ssgsea/xiaoman/expression/annot/%s.annot",dataset.prefix)

output.dir = sprintf("~/project/deeplearning/icb/data/%s", dataset.prefix)

pca_obj.RData = "/homes6/asahu/project/deeplearning/icb/data/tcga.blca/neoantigen/pca_obj.RData"

fix_patient_name =F; ICB_dataset =T

# precog dataset 
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Avin/Precog_GSVA.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_Precog_oxphos.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm"

output.dir = "~/project/deeplearning/icb/data/tcga.blca/tcga.blca.pca/"

#######
######
#######
######
#######
######
#######
######
data_dir="/homes6/asahu/project/deeplearning/icb/data/ssgsea/datasets/"
model_dir="/homes6/asahu/project/deeplearning/icb/data/ssgsea/models/"




# python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation/datasets
python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/datasets --prefix "tcga"
python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation15Feb/datasets_list.txt --model_dir  ~/Dropbox/project/code/deeplearning/icb/results/simulation15Feb/

python train.py --data_dir $data_dir --model_dir $model_dir

python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/datasets_list.txt --prefix "tcga"
# simulation from pixie 

python train.py --data_dir ../results/simulation5Feb/datasets_list.txt --model_dir ../results/simulation5Feb/


python train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43

ipython -pdb train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43


python train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_embedding8


python train.py --prefix tcga --data_dir experiments/tcga_feb9/datasets_list.txt  --model_dir experiments/tcga_feb9

python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/tcga.oxphos/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/tcga.oxphos/


python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/ --restore_file tcga_survial

python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/survival/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/survival/ --restore_file tcga_survial

python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_icb/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_icb/ 

python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_icb/icb/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_icb/icb/ --restore_file tcga_icb_last 


python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_only/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_only/

python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_only/icb_survival/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/tcga_only/icb_survival/ --restore_file intialize

#precog
python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/
#tcga
python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/tcga/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/tcga/
python  train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/combined/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/combined/

python train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/tcga.blca/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/tcga.blca 

python train.py --data_dir ../data/tcga.blca/genentech.pca/datasets_list.txt --model_dir ../data/tcga.blca/genentech.pca --restore_file ../data/tcga.blca/pca/tensorboardLog/20190307-140948/best.pth.tar



python evaluate.py  --data_dir /Users/avi/Dropbox/project/code/deeplearning/icb/data/tcga.blca/genentech.pca/eval/datasets_list.txt  --model_dir /Users/avi/Dropbox/project/code/deeplearning/icb/data/tcga.blca/pca/tensorboardLog/20190307-140948/ --restore_file /Users/avi/Dropbox/project/code/deeplearning/icb/data/tcga.blca/pca/tensorboardLog/20190307-140948/best.pth.tar


python train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/tcga/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/tcga/ --restore_file ~/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/tcga/last.best.pthr.tar



python train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/pretrain_tcga/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/pretrain_tcga/ 
--restore_file ~/project/deeplearning/icb/data/pancancer_all_immune/precog.oxphos/norm/tcga/my.best.pth.tar



python train.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/tcga.oxphos/tcga_genentech/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/tcga.oxphos/tcga_genentech/

python evaluate.py  --data_dir  /homes6/asahu/project/deeplearning/icb/data/tcga.oxphos/tcga_genentech/datasets_list.txt --model_dir /homes6/asahu/project/deeplearning/icb/data/tcga.oxphos/tcga_genentech/ --restore_file best


python train.py  --data_dir  ../data/genentech.tpm/datasets_list.txt --model_dir ../data/genentech.tpm/. --tensorboard_prefix only_auc_ 


 python  train.py  --data_dir  ../data/genentech.tpm/genentech.pca.tpm.phenotypes.32/datasets_list.txt --model_dir ../data/genentech.tpm/genentech.pca.tpm.phenotypes/ --tensorboard_prefix "32_"


  --restore_file ../data/genentech.tpm/genentech.pca.phenotypes/tensorboardLog/surv_20190310-214230/best.pth.tar

## Neoantigen
 
python train.py  --data_dir  ../data/tcga/PCA/datasets_list.txt --model_dir ../data/tcga/PCA/. 

python  evaluate.py  --data_dir  ../data/tcga/PCA/datasets_list.txt --model_dir ../data/tcga/PCA/.  --restore_file  ../data/tcga/PCA/./tensorboardLog/20190322-014624/best.pth.tar

python  evaluate.py  --data_dir   ../data/tcga/PCA/prediction_genentech/datasets_list.txt --model_dir ../data/tcga/PCA/prediction_genentech/.  --restore_file  ../data/tcga/PCA/./tensorboardLog/20190322-014624/best.pth.tar


#Neoantigen in BCB 
python train.py  --data_dir  ../data/tcga.blca/neoantigen/datasets_list.txt --model_dir ../data/tcga.blca/neoantigen/. 

python evaluate.py  --data_dir  ../data/tcga/PCA/datasets_list.txt --model_dir ../data/tcga/PCA/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar

python  evaluate.py  --data_dir  ../data/genentech.tpm/Neoantigen/evaluate/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar

python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar


python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar --tensorboard_prefix pretrain_tcga_



## Response prediction
python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen.imputed/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen.imputed/. --tensorboard_prefix auc_

python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen.imputed/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen.imputed/. --tensorboard_prefix auc_neo_ml_ --hyper_param "input_indices=range(50)"


python  -m ipdb train.py  --data_dir  ../data/genentech.tpm/Neoantigen.imputed/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen.imputed/. --tensorboard_prefix outputs_10_ --hyper_param "params.binary_phenotype_indices=[45,1,44,47,46,41];params.continuous_phenotype_indices=[0,42,43,49,50]"

python  -m ipdb train.py  --data_dir  ../data/genentech.tpm/Neoantigen.imputed/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen.imputed/. --tensorboard_prefix outputs_10_ --hyper_param "params.continuous_phenotype_indices=[0,42,43,49,50]"
# run run_hyper.R  with ii = 57

## mGC_PD1_Kim
python  evaluate.py \
 --data_dir  ../data/mGC_PD1_Kim/datasets_list.txt \
 --model_dir ../data/mGC_PD1_Kim/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.binary_phenotype_indices=[42,43,45,1,44,47,46,41];params.continuous_phenotype_indices=[0,50,49,55,54,56,53,17,51,52,29,13,14,31,10];params.input_indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,60,56]"


# imputing neoantigen etc. from tcga model 
python  evaluate.py \
 --data_dir  ../data/mGC_PD1_Kim/Neoantigen/datasets_list.txt \
 --model_dir ../data/mGC_PD1_Kim/Neoantigen/. \
 --restore_file ../data/tcga/PCA/tensorboardLog/20190322-015246/best.pth.tar


python  evaluate.py \
 --data_dir  ../data/mGC_PD1_Kim/imputed/datasets_list.txt \
 --model_dir ../data/mGC_PD1_Kim/imputed/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.binary_phenotype_indices=[42,43,45,1,44,47,46,41];params.continuous_phenotype_indices=[0,50,49,55,54,56,53,17,51,52,29,13,14,31,10];params.input_indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,60,56]"





## RCC_PD1_Miao
python  evaluate.py \
 --data_dir  ../data/RCC_PD1_Miao/datasets_list.txt \
 --model_dir ../data/RCC_PD1_Miao/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar

 # Mel_CTLA4_VanAllen 


 python  evaluate.py \
 --data_dir  ../data/Mel_CTLA4_VanAllen/datasets_list.txt \
 --model_dir ../data/Mel_CTLA4_VanAllen/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.cuda=False" 

python  evaluate.py \
 --data_dir  ../data/Mel_CTLA4_VanAllen/Neoantigen/datasets_list.txt \
 --model_dir ../data/Mel_CTLA4_VanAllen/Neoantigen/. \
 --restore_file ../data/tcga/PCA/tensorboardLog/20190322-015246/best.pth.tar

 python  evaluate.py \
 --data_dir  ../data/Mel_CTLA4_VanAllen/imputed/datasets_list.txt \
 --model_dir ../data/Mel_CTLA4_VanAllen/imputed/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.binary_phenotype_indices=[42,43,45,1,44,47,46,41];params.continuous_phenotype_indices=[0,50,49,55,54,56,53,17,51,52,29,13,14,31,10];params.input_indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,60,56]"




 # Mel_CTLA4_VanAllen 


 python  evaluate.py \
 --data_dir  ../data/Mel_PD1_Hugo/datasets_list.txt \
 --model_dir ../data/Mel_PD1_Hugo/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.cuda=False" 


python  evaluate.py \
 --data_dir  ../data/Mel_PD1_Hugo/Neoantigen/datasets_list.txt \
 --model_dir ../data/Mel_PD1_Hugo/Neoantigen/. \
 --restore_file ../data/tcga/PCA/tensorboardLog/20190322-015246/best.pth.tar



## cancer type specific prediction 
## Note TCGA data is being used in the code 
## 



# imputing neoantigen etc. from tcga model 
python  evaluate.py \
 --data_dir  ../data/tcga/PCA/Neoantigen/datasets_list.txt \
 --model_dir ../data/tcga/PCA/Neoantigen/. \
 --restore_file ../data/tcga/PCA/tensorboardLog/20190322-015246/best.pth.tar


python  evaluate.py \
 --data_dir  ../data/tcga/PCA/imputed/datasets_list.txt \
 --model_dir ../data/tcga/PCA/imputed/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.binary_phenotype_indices=[42,43,45,1,44,47,46,41];params.continuous_phenotype_indices=[0,50,49,55,54,56,53,17,51,52,29,13,14,31,10];params.input_indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,60,56]"


## genetech response prediction
python  evaluate.py \
 --data_dir  ../data/genentech.tpm/Neoantigen.imputed/evaluate/datasets_list.txt \
 --model_dir ../data/genentech.tpm/Neoantigen.imputed/evaluate/. \
 --restore_file ~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/tensorboardLog/output_22/inps_61_5720190321-104715/best.pth.tar \
 --hyper_param "params.binary_phenotype_indices=[42,43,45,1,44,47,46,41];params.continuous_phenotype_indices=[0,50,49,55,54,56,53,17,51,52,29,13,14,31,10];params.input_indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,60,56]"



 ## version > 0.3
 # ../data/tcga/neoantigen.v2

python train.py  --data_dir  ../data/tcga/neoantigen.v2/datasets_list.txt --model_dir ../data/tcga/neoantigen.v2/. --tensorboard_prefix msi_only_ 

python train.py  --data_dir  ../data/ya/datasets_list.txt --model_dir ../data/ya/.

# attention model 


python train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/.


python evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/.  --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/20190407-213630/best.pth.tar



python evaluate_attention.py  --data_dir  ../data/genentech.tpm/neoantigen.v2/impute/datasets_list.txt --model_dir ../data/genentech.tpm/neoantigen.v2/impute/.  --restore_file  ../data/tcga/neoantigen.v2/attention/tensorboardLog/20190412-151621/best.pth.tar


python train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/.


python train.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/.

## ya's model 
python train.py  --data_dir  ../data/ya/datasets_tsne_list.txt --model_dir ../data/ya/.

## imputation 
python evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_tsne_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/.  --restore_file  ../data/tcga/neoantigen.v2/attention/tensorboardLog/20190422-080613/best.pth.tar --type_file="['']"

python evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/survival20190422-081755/datasets_impute_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/survival20190422-081755/.  --restore_file  ../data/tcga/neoantigen.v2/attention/tensorboardLog/survival20190422-081755/best.pth.tar --type_file="all"


python train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/genentech.imputed/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/genentech.imputed/.

python train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/. --restore_file ~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/tensorboardLog/SGD_20190426-131226//best.pth.tar

python train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/. --restore_file ~/Dropbox/project/code/deeplearning/icb/data/tcga/neoantigen.v2/attention/tcga.imputed/tensorboardLog/SGD_20190426-131226//best.pth.tar

python evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/. --restore_file ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/best.pth.tar

python evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/datasets_genentech_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/. --restore_file ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/best.pth.tar




# kindney

python train.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/. --tensorboard_prefix  new_residual_

python -m ipdb  evaluate.py   --data_dir ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/datasets_tsne_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/best.pth.tar

python -m ipdb evaluate.py   --data_dir ../data/RCC_PD1_Miao/datasets_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/best.pth.tar


python evaluate.py   --data_dir ../data/RCC_PD1_Miao/datasets_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/best.pth.tar



python  evaluate.py   --data_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/datasets_genentech_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190515-114038/best.pth.tar

## evaluating kidney cancer 
python train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/drug/.


python  train_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/.

python  evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/. --restore_file  ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/tensorboardLog/three_var_20190622-191219/best.pth.tar

python  evaluate_attention.py  --data_dir  ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/. --restore_file  ../data/tcga/neoantigen.v2/attention/tcga.imputed/brca/tensorboardLog/20190623-230228/best.pth.tar


# scrRNA

python train.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/. 

python train.py  --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/. --tensorboard_prefix  restored_ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190625-233202/best.pth.tar

python evaluate.py   --data_dir  ../data/tcga/neoantigen.v2/attention/datasets_val_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190625-041025/best.pth.tar


python evaluate.py   --data_dir ../data/Getz_scRNA/datasets_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190625-041025/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190625-041025/best.pth.tar

# in pixie 

python evaluate.py   --data_dir ../data/Getz_scRNA/datasets_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190625-233202/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/new_residual_20190625-233202/best.pth.tar

python evaluate.py   --data_dir ../data/Getz_scRNA/datasets_list.txt --model_dir  ../data/tcga/neoantigen.v2/attention/tensorboardLog/restored_20190627-185716/ --restore_file ../data/tcga/neoantigen.v2/attention/tensorboardLog/restored_20190627-185716/best.pth.tar


dca /homes6/asahu/project/deeplearning/icb/data/Getz_scRNA/data/count.mat.txt /homes6/asahu/project/deeplearning/icb/data/Getz_scRNA/data/dca_impute/


# tar -xf /rna-seq/dbGAP_tar/Post_P16_CD45+.tar --wildcards --no-anchored '*.tab'  -C /rna-seq/counts/.  

# scRNA version deepImmune 

python train.py  --data_dir  ../data/tcga//scrna.v1/datasets_tsne_list.txt --model_dir ../data/tcga//scrna.v1/. 

python evaluate.py  --data_dir  ../data/Getz_scRNA/dca/datasets_list.txt --model_dir ../data/Getz_scRNA/dca/.  --restore_file ../data/tcga/scrna.v1/saved_model/20190715-115243/best.pth.tar


python evaluate.py  --data_dir  ../data/tcga//scrna.v1/datasets_test_list.txt --model_dir ../data/tcga/scrna.v1/saved_model/20190715-115243/ --restore_file ../data/tcga/scrna.v1/saved_model/20190715-115243/best.pth.tar

python train.py  --data_dir  ../data/tcga//scrna.v1/datasets_tsne_list.txt --model_dir ../data/tcga//scrna.v1/.  --restore_file ../data/tcga/scrna.v1/tensorboardLog/20190715-115243/best.pth.tar --tensorboard_prefix boosting_

python train.py  --data_dir  ../data/tcga//scrna.v2/datasets_tsne_list.txt --model_dir ../data/tcga//scrna.v2/. --restore_file ../data/tcga/scrna.v2/tensorboardLog/20190717-195353/best.pth.tar --tensorboard_prefix boosting_

python train.py  --data_dir  ../data/tcga//scrna.v2/datasets_tsne_list.txt --model_dir ../data/tcga//scrna.v2/. --tensorboard_prefix indel_only_1ouptput_

python train.py  --data_dir  ../data/tcga//scrna.v3/datasets_tsne_list.txt --model_dir ../data/tcga//scrna.v3/. 
python train.py  --data_dir  ../data/tcga/scrna.v4.pcs/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.pcs/. 


python train.py  --data_dir  ../data/tcga/scrna.v4.genes/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.genes/.  --tensorboard_prefix no_pipeline_


python train.py  --data_dir  ../data/tcga/scrna.v4.allgenes.nopcs/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.allgenes.nopcs/.  --tensorboard_prefix no_pipeline_



python evaluate.py  --data_dir  ../data/Getz_scRNA/scrna.v4.genes/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/best.pth.tar

python evaluate.py  --data_dir  ../data/tcga/scrna.v4.genes/datasets_test_list.txt --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/best.pth.tar

CUDA_VISIBLE_DEVICES=2 python evaluate.py  --data_dir  ../data/Getz_scRNA/scrna.v4.genes/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/epoch-52.pth.tar

CUDA_VISIBLE_DEVICES=2 python evaluate.py  --data_dir  ../data/tcga/scrna.v4.genes/datasets_test_list.txt  --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/epoch-52.pth.tar


CUDA_VISIBLE_DEVICES=2 python evaluate.py  --data_dir  ../data/oxphos/scrna.v4.genes/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/no_pipeline_20190724-201216/epoch-52.pth.tar





### TCR Amino acid training 

CUDA_VISIBLE_DEVICES=1 python train.py  --data_dir  ../data/Getz_scRNA/TCR.AA/datasets_tsne_list.txt --model_dir ../data/Getz_scRNA/TCR.AA/. --tensorboard_prefix emb_n_immuneFactors_ 

python evaluate.py  --data_dir  ../data/Getz_scRNA/TCR.AA/datasets_test_list.txt --model_dir ../data/Getz_scRNA/TCR.AA/autoencoder_20190820-140204/. --restore_file ../data/Getz_scRNA/TCR.AA/tensorboardLog/autoencoder_20190820-140204/best.pth.tar


# with  VAE 
# mac 
python train.py  --data_dir  ../data/tcga/neoantigen.v2/attention//datasets_tsne_list.txt --model_dir ../data/tcga/neoantigen.v2/attention/


python train.py  --data_dir  ../data/tcga/scrna.v4.genes/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.genes/.  --tensorboard_prefix nopipeline_vae_


CUDA_VISIBLE_DEVICES=2 python evaluate.py  --data_dir  ../data/Getz_scRNA/scrna.v4.genes/datasets_tsne_list.txt --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/nopipeline_vae_20190819-161104/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/nopipeline_vae_20190819-161104/epoch-142.pth.tar

CUDA_VISIBLE_DEVICES=2 python evaluate.py  --data_dir  ../data/tcga/scrna.v4.genes/datasets_test_list.txt --model_dir ../data/tcga/scrna.v4.genes/tensorboardLog/nopipeline_vae_20190819-161104/. --restore_file  ../data/tcga/scrna.v4.genes/tensorboardLog/nopipeline_vae_20190819-161104/epoch-142.pth.tar

