data_dir="/homes6/asahu/project/deeplearning/icb/data/ssgsea/datasets/"
model_dir="/homes6/asahu/project/deeplearning/icb/data/ssgsea/models/"

# python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation/datasets
python train.py --data_dir ~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/datasets
python train.py --data_dir $data_dir --model_dir $model_dir


python train.py --prefix tcga --data_dir ~/local_data/processed/datasets --model_dir /homes6/asahu/project/deeplearning/icb/deepImmune/experiments/tcga_43