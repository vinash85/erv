###
#######
#processing of datasets
######
# ICB dataset 
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Avin/ICB_GSVA.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_ICB.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/all_icb"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech.imputed"


output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/genetech.pca/"
fix_patient_name =F; ICB_dataset =T

# precog dataset 
dataset_ssgsea = "/liulab/asahu/data/ssgsea/xiaoman/Avin/Precog_GSVA.txt"
pathway_order = "/liulab/asahu/data/ssgsea/xiaoman/ssgsea.order_tcga.txt"
dataset_phenotype = "/liulab/asahu/data/ssgsea/xiaoman/Avin/clinical_Precog.txt"
phenotype_order = "/liulab/asahu/data/ssgsea/xiaoman/processed/tcga_phenotypes.RData"
output.dir = "~/project/deeplearning/icb/data/pancancer_all_immune/precog"


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

/homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43/*.tar

python train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_embedding8


python train.py --prefix tcga --data_dir experiments/tcga_feb9/datasets_list.txt  --model_dir experiments/tcga_feb9