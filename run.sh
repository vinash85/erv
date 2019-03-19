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
output.dir = "~/project/deeplearning/icb/data/tcga.brca/PCA/"

output.dir = "~/project/deeplearning/icb/data/tcga.blca/neoantigen/"

# ICB dataset

dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Avin/ICB_GSVA.txt"
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Genetech_expression_TPM.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_ICB.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/all_icb"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed.same.survival/pca/"
output.dir = "~/project/deeplearning/icb/data/genentech.tpm/Neoantigen/"
output.dir = "~/project/deeplearning/icb/data/genentech.tpm/Neoantigen.imputed/"


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


python train.py  --data_dir  ../data/genentech.tpm/genentech.pca.tpm.phenotypes/datasets_list.txt --model_dir ../data/genentech.tpm/genentech.pca.tpm.phenotypes/. 

--tensorboard_prefix ""

 python  train.py  --data_dir  ../data/genentech.tpm/genentech.pca.tpm.phenotypes.32/datasets_list.txt --model_dir ../data/genentech.tpm/genentech.pca.tpm.phenotypes/ --tensorboard_prefix "32_"


  --restore_file ../data/genentech.tpm/genentech.pca.phenotypes/tensorboardLog/surv_20190310-214230/best.pth.tar

## Neoantigen
 
python train.py  --data_dir  ../data/tcga/PCA/datasets_list.txt --model_dir ../data/tcga/PCA/.
python evaluate.py  --data_dir  ../data/tcga/PCA/datasets_list.txt --model_dir ../data/tcga/PCA/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar
python  evaluate.py  --data_dir  ../data/genentech.tpm/Neoantigen/evaluate/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar
python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar


python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen/. --restore_file ../data/tcga/PCA/tensorboardLog/20190315-150358/best.pth.tar --tensorboard_prefix pretrain_tcga_


python  train.py  --data_dir  ../data/genentech.tpm/Neoantigen.imputed/datasets_list.txt --model_dir ../data/genentech.tpm/Neoantigen.imputed/. --tensorboard_prefix auc_

