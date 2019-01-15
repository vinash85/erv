# read the embedding model and variable
import torch
import numpy as np

# model_dir = "./experiments/base_model/"
model_dir = "./experiments/tcga_43/"

checkpoint = model_dir + "best.pth.tar"
checkpoint = torch.load(checkpoint)
# model.load_state_dict(checkpoint['state_dict'])
embedding_params = checkpoint['embedding_state_dict']
inp_embedding_weight = embedding_params['fc1.weight']
# note weight matrix AxB output_dim x input_dim

outputs_params = checkpoint['outputs']
linear_output_weights = outputs_params['linear1.weight']
binary_output_weights = outputs_params['linear2.weight']

np.savetxt(model_dir + "linear_output_weights.txt", linear_output_weights, delimiter=",")
