data_dir="/homes6/asahu/project/deeplearning/icb/data/ssgsea/datasets/"
model_dir="/homes6/asahu/project/deeplearning/icb/data/ssgsea/models/"

# python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation/datasets
python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/datasets --prefix "tcga"
python train.py --data_dir $data_dir --model_dir $model_dir

python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/datasets_list.txt --prefix "tcga"

python train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43

ipython -pdb train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43

/homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43/*.tar

python train.py --prefix tcga --data_dir ~/local_data/processed/datasets/ --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_embedding8