from lifelines.utils import concordance_index
import numpy as np


def concordance_metric(predicted_risk, survival):
    # calculate the concordance index
    survival_time, censor = survival[:, 0], survival[:, 1]
    epsilon = 0.001
    partial_hazard = np.exp(-(predicted_risk + epsilon))
    censor = censor.astype(int)
    ci = concordance_index(survival_time, partial_hazard, censor)
    return ci
