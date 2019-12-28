import json
import logging
import os
import shutil
import torch
import jstyleson
import numpy as np
import math
from past.builtins import basestring
import umap
from sklearn.datasets import fetch_openml
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


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
        for key, value in params.items():
            if isinstance(value, basestring):
                if value[0] == "e":
                    params[key] = eval(value[1:])
        return params

    @property
    def dict(self):
        """Gives dict-like access to Params instance by `params.dict['learning_rate']"""
        return self.__dict__


class RunningAverage():
    """A simple class that maintains the running average of a quantity

    Example:
    ```
    loss_avg = RunningAverage()
    loss_avg.update(2)
    loss_avg.update(4)
    loss_avg() = 3
    ```
    """

    def __init__(self):
        self.steps = 0
        self.total = 0

    def update(self, val):
        self.total += val
        self.steps += 1

    def __call__(self):
        return self.total / float(self.steps)


def set_logger(log_path):
    """Set the logger to log info in terminal and file `log_path`.

    In general, it is useful to have a logger so that every output to the terminal is saved
    in a permanent file. Here we save it to `model_dir/train.log`.

    Example:
    ```
    logging.info("Starting training...")
    ```

    Args:
        log_path: (string) where to log
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        # Logging to a file
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(logging.Formatter('%(asctime)s:%(levelname)s: %(message)s'))
        logger.addHandler(file_handler)

        # Logging to console
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(stream_handler)


def save_dict_to_json(d, json_path):
    """Saves dict of floats in json file

    Args:
        d: (dict) of float-castable values (np.float, int, float, etc.)
        json_path: (string) path to json file
    """
    with open(json_path, 'w') as f:
        # We need to convert the values to float for json (it doesn't accept np.array, np.float, )
        try:
            d = {k: float(v) for k, v in d.items()}
        except:
            d = d

        jstyleson.dump(d, f, indent=4)


def save_checkpoint(state, is_best, checkpoint, epoch=None):
    """Saves model and training parameters at checkpoint + 'last.pth.tar'. If is_best==True, also saves
    checkpoint + 'best.pth.tar'

    Args:
        state: (dict) contains model's state_dict, may contain other keys such as epoch, optimizer state_dict
        is_best: (bool) True if it is the best model seen till now
        checkpoint: (string) folder where parameters are to be saved
    """
    filepath = os.path.join(checkpoint, 'last.pth.tar')
    if not os.path.exists(checkpoint):
        print("Checkpoint Directory does not exist! Making directory {}".format(checkpoint))
        os.mkdir(checkpoint)
    else:
        print("Checkpoint Directory exists! ")
    torch.save(state, filepath)
    if epoch is not None:
        shutil.copyfile(filepath, os.path.join(checkpoint, 'epoch-{}.pth.tar'.format(epoch)))

    if is_best:
        shutil.copyfile(filepath, os.path.join(checkpoint, 'best.pth.tar'))


def load_checkpoint(checkpoint, embedding_model, outputs, embedding_optimizer=None, outputs_optimizer=None):
    """Loads model parameters (state_dict) from file_path. If optimizer is provided, loads state_dict of
    optimizer assuming it is present in checkpoint.

    Args:
        checkpoint: (string) filename which needs to be loaded
        model: (torch.nn.Module) model for which the parameters are loaded
        optimizer: (torch.optim) optional: resume optimizer from checkpoint
    """
    if not os.path.exists(checkpoint):
        raise "File doesn't exist {}"
    checkpoint = torch.load(checkpoint)
    embedding_model.load_state_dict(checkpoint['embedding_state_dict'])
    outputs.load_state_dict(checkpoint['outputs_state_dict'])

    if embedding_optimizer:
        embedding_optimizer.load_state_dict(checkpoint['embedding_optim_dict'])
    if outputs_optimizer:
        outputs_optimizer.load_state_dict(checkpoint['outputs_optim_dict'])

    return checkpoint


def load_checkpoint_attn(checkpoint, models, params, optimizers=None):
    """Loads model parameters (state_dict) from file_path. If optimizer is provided, loads state_dict of
    optimizer assuming it is present in checkpoint.

    Args:
        checkpoint: (string) filename which needs to be loaded
        model: (torch.nn.Module) model for which the parameters are loaded
        optimizer: (torch.optim) optional: resume optimizer from checkpoint
    """
    if not os.path.exists(checkpoint):
        raise "File doesn't exist {}"

    if params.cuda:
        checkpoint = torch.load(checkpoint)
    else:
        checkpoint = torch.load(checkpoint, map_location='cpu')

    models[0].load_state_dict(checkpoint['embedding_state_dict'])
    models[1].load_state_dict(checkpoint['attention_state_dict'])
    models[2].load_state_dict(checkpoint['outputs_state_dict'])

    if optimizers:
        optimizers[0].load_state_dict(checkpoint['embedding_optim_dict'])
        optimizers[1].load_state_dict(checkpoint['attention_optim_dict'])
        optimizers[2].load_state_dict(checkpoint['outputs_optim_dict'])

    return checkpoint


def reshape_toimage(mat, dim1=None):
    size1 = mat.size
    if dim1 is None:
        dim1 = math.ceil(size1**0.5 + 0.5)
    dim2 = math.ceil(size1 / dim1 + 0.5)
    mat = np.pad(mat, (0, dim1 * dim2 - size1), "constant", constant_values=(0, 0))
    mat = 1 - mat.reshape(1, dim1, dim2)
    return mat


def create_save_umap(mat, color, epoch, savedir=None, title=None):
    '''
    Create umap and save figure
    '''

    embedding = umap.UMAP(n_neighbors=15,
                          min_dist=0.1,
                          metric='correlation').fit_transform(mat)
    # col = color.astype(int)
    df = pd.DataFrame(dict(color=color))
    fig = plt.figure()
    # plt.scatter(embedding[:, 0], embedding[:, 1])
    # sns.lmplot('dim1', 'dim2', hue='color', data=df, fit_reg=False)
    # Unique category labels: 'D', 'F', 'G', ...
    color_labels = df['color'].unique()

    # List of RGB triplets
    rgb_values = sns.color_palette("Set2", len(color_labels))

    # Map label to RGB
    color_map = dict(zip(color_labels, rgb_values))

    # Finally use the mapped values
    plt.scatter(embedding[:, 0], embedding[:, 1], c=df['color'].map(color_map))
    if title is not None:
        plt.title(title)
    if savedir is not None:
        plt.savefig('{}/epoch-{}.png'.format(savedir, epoch), format='png', dpi=100)
    return fig
