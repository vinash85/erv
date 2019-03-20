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
