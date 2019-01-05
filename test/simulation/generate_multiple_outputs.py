import numpy as np
# import matplotlib.pyplot as plt
import sys
import pandas as pd
import random 
sys.path.insert(0, '~/Dropbox/project/code/deeplearning/antigen_recognition/src')
sys.path.append('/Users/avi/Dropbox/project/code/deeplearning/antigen_recognition/src')
import r2python
numfeat = 100
numSample = 2000

feature = np.random.normal(0, 1, numfeat * numSample)
featureMat = np.mat(feature)
featureMat = featureMat.reshape(numSample, numfeat)

risk1 = .3 * featureMat[:, 90] - .7 * featureMat[:, 94]
risk = np.squeeze(np.asarray(risk1))
risk_noise = risk + np.random.normal(0, np.std(risk) / 2, numSample)
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

survival_final.to_csv("~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/survival.train.txt", sep='\t', index=False)
# numpy.savetxt("../results/simulation/survival.txt", a, delimiter=",")
feature_df = pd.DataFrame(featureMat)

def create_outputs(featureMat, linear=True):

	# choose two variable
	out = featureMat[:,random.randint(0, 99)] * random.uniform(-1,1) +  featureMat[:,random.randint(0, 99)] * random.uniform(-1,1) 
	# out1 = out 
	if not linear:
		for i in range(len(out)): 
			if out[i] > 0 : 
				out[i] = 1
			else :
				out[i] = 0
	return out

linear_output_size = 10
binary_output_size = 10

# outputs = np.zeros( numSample , linear_output_size + binary_output_size + 1)

outputs = survival_final

for i in range(linear_output_size-1) :
 	tag = "linear" + str(i)
	outputs[tag] = create_outputs(featureMat)

for i in range(binary_output_size) :
 	tag = "binary" + str(i)
	outputs[tag] = create_outputs(featureMat)


outputs.to_csv("~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/outputs.train.txt", sep='\t', index=False)
feature_df.to_csv("~/Dropbox/project/code/deeplearning/icb/results/simulation5jan/survival.data.txt", sep='\t', index=False)
