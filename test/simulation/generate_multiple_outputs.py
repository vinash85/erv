import numpy as np
# import matplotlib.pyplot as plt
import sys
import pandas as pd
import random
import math
import os
# sys.path.insert(0, '~/Dropbox/project/code/deeplearning/antigen_recognition/src')
# sys.path.append('/Users/avi/Dropbox/project/code/deeplearning/antigen_recognition/src')
# sys.path.append('/Users/avi/Dropbox/project/code/deeplearning/antigen_recognition/src')
import r2python
# dataset_dir = "~/Dropbox/project/code/deeplearning/icb/results/simulation5Feb/"
dataset_dir = "../results/simulation15Feb/"
os.system("mkdir " + dataset_dir)
numfeat = 100
numSample = 20000
noise_level = 2.0  # amount of noise

feature = np.random.normal(0, 1, numfeat * numSample)
featureMat = np.mat(feature)
featureMat = featureMat.reshape(numSample, numfeat)

risk1 = sum([featureMat[:, random.randint(0, numfeat - 1)] * random.uniform(-1, 1) for ii in range(5)])
risk = np.squeeze(np.asarray(risk1))
risk_noise = risk + np.random.normal(0, np.std(risk) * noise_level, numSample)
risk_noise = risk_noise - min(risk_noise)
survival = risk_noise
censor = np.random.choice([0, 1], numSample)


def create_survival(cc, ss):
    out = ss
    if(cc == 0):
        out = np.random.uniform(0, 1) * ss
    return(out)


survival_censored = [create_survival(cc, ss) for cc, ss in zip(censor, survival)]
# np.corrcoef(np.asarray(survival_censored), np.squeeze(np.asarray(featureMat[:, 90])))
aa = np.asarray([survival_censored, censor])
aa = aa.transpose([1, 0])
survival_final = pd.DataFrame(data=aa, columns=["times", "status"])

survival_final.to_csv(dataset_dir + "survival.train.txt", sep='\t', index=False)
# numpy.savetxt("../results/simulation/survival.txt", a, delimiter=",")
feature_df = pd.DataFrame(featureMat)


def create_outputs(featureMat, linear=True, add_nan=False):

        # choose two variable
    out1 = featureMat[:, random.randint(0, 99)] * random.uniform(-1, 1) + featureMat[:, random.randint(0, 99)] * random.uniform(-1, 1)
    noise = np.random.normal(0, np.std(out1) * noise_level, numSample)
    out = out1 + np.reshape(noise, [numSample, 1])
    if not linear:
        for i in range(len(out)):
            if out[i] > 0:
                out[i] = 1
            else:
                out[i] = 0
    # add nan to 5% of data
    if add_nan:
        arr = np.arange(len(out))
        np.random.shuffle(arr)
        sel = arr[:int(math.ceil(len(out) * .05))]
        out[sel] = np.nan
    return out

linear_output_size = 2
binary_output_size = 1

# outputs = np.zeros( numSample , linear_output_size + binary_output_size + 1)

outputs = survival_final

for i in range(linear_output_size - 1):
    tag = "linear" + str(i)
    outputs[tag] = create_outputs(featureMat, add_nan=True)

for i in range(binary_output_size):
    tag = "binary" + str(i)
    outputs[tag] = create_outputs(featureMat, linear=False)


outputs.to_csv(dataset_dir + "outputs.train.txt", sep='\t', index=False, na_rep='NaN')
feature_df.to_csv(dataset_dir + "features.txt", sep='\t', index=False)


# split 7 2 1
dataset_dir1 = dataset_dir + "datasets/"
os.system("mkdir " + dataset_dir1)

feature_df.iloc[:int(.7 * len(outputs)), :].to_csv(dataset_dir1 + "ssgsea_train.txt", sep='\t', index=False, na_rep='NaN')
outputs.iloc[:int(.7 * len(outputs)), :].to_csv(dataset_dir1 + "phenotype_train.txt", sep='\t', index=False, na_rep='NaN')

feature_df.iloc[int(.7 * len(outputs)): int(.9 * len(outputs)), :].to_csv(dataset_dir1 + "ssgsea_val.txt", sep='\t', index=False, na_rep='NaN')
outputs.iloc[int(.7 * len(outputs)): int(.9 * len(outputs)), :].to_csv(dataset_dir1 + "phenotype_val.txt", sep='\t', index=False, na_rep='NaN')

feature_df.iloc[int(.9 * len(outputs)):, :].to_csv(dataset_dir1 + "ssgsea_test.txt", sep='\t', index=False, na_rep='NaN')
outputs.iloc[int(.9 * len(outputs)):, :].to_csv(dataset_dir1 + "phenotype_test.txt", sep='\t', index=False, na_rep='NaN')
