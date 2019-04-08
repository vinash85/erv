# debuging that loss is zero for one of the encoder ouput the derivative is also zero
params.loss_excluded_from_training
list(outputs.parameters())
for p in outputs.named_parameters():
    print(p)

cc = outputs.parameters().next()


import json
import logging
import os
import shutil
import torch
import jstyleson


class Params():
    """Class that loads hyperparameters from a json file.

    Example:
    ```
    params = Params(json_path)
    print(params.learning_rate)
    params.learning_rate = 0.5  # change the value of learning_rate in params
    ```
    """

    def __init__(self, json_path):
        with open(json_path) as f:
            params = jstyleson.load(f)
            params = self.eval_string(params)
            self.__dict__.update(params)

    def save(self, json_path):
        with open(json_path, 'w') as f:
            jstyleson.dump(self.__dict__, f, indent=4)

    def update(self, json_path):
        """Loads parameters from json file"""
        with open(json_path) as f:
            params = jstyleson.load(f)
            self.__dict__.update(params)

    def eval_string(self, params):
        for key, value in params.iteritems():
            if isinstance(value, basestring):
                if value[0] == "e":
                    params[key] = eval(value[1:])
        return params

    @property
    def dict(self):
        """Gives dict-like access to Params instance by `params.dict['learning_rate']"""
        return self.__dict__


m = nn.Sigmoid()
loss = nn.BCELoss()
input = torch.randn(3, requires_grad=True)
target = torch.empty(3).random_(2)
output = loss(m(input), target)
output.backward()


# for debugging tensor for nans

import torch
from torch.autograd.variable import Variable
import torch.nn.functional as F

a = Variable(torch.FloatTensor(16, 3, 64, 64).zero_(), requires_grad=True)

# b = a.std(dim=0).mean()  # NaN gradients
b = (a - a.mean(dim=0, keepdim=True)).norm(p=2, dim=0).mean() / (15**0.5)  # 0 gradients

b.backward()

print(a.grad.data)


if not all([isfinite(p.grad).all() for p in outputs.parameters()]):
    tracer()

# get intermediate ouputs
# https://forums.fast.ai/t/pytorch-best-way-to-get-at-intermediate-layers-in-vgg-and-resnet/5707
# res50_conv = nn.Sequential(*list(res50_model.children())[:-2])
yy = [xx for xx in embedding_model.named_parameters()]
aa = nn.Sequential(*list(embedding_model.children())[:2])
