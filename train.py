"""Train the model"""


import torch
# import torch.nn as nn
# import torchvision
# import torchvision.transforms as transforms
# from ... import data_generator as gn
# import data_generator_pytorch as gn
# import datetime
# import time


import argparse
import logging
import os
import numpy as np
import torch
import torch.optim as optim
from torch.autograd import Variable
from tqdm import tqdm

import utils
import model.net as net
# import model.data_loader as data_loader
import model.data_generator as data_generator
from evaluate import evaluate

parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', default='data/64x64_SIGNS', help="Directory containing the dataset")
parser.add_argument('--model_dir', default='experiments/base_model', help="Directory containing params.json")
parser.add_argument('--restore_file', default=None,
                    help="Optional, name of the file in --model_dir containing weights to reload before \
                    training")  # 'best' or 'train'

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
# cudnn.benchmark = True

# # Hyper parameters
# num_epochs = 200
# num_classes = 1
# batch_size = 100
# learning_rate = 0.001
# hidden_size = 256
# now = datetime.datetime.now
# t = now()
# timestr = time.strftime("%Y%m%d_%H%M")

# train_batch_size = 256


def train(model, optimizer, loss_fn, dataloader, metrics, params):
    """Train the model on `num_steps` batches
    Args:
        model: (torch.nn.Module) the neural network
        optimizer: (torch.optim) optimizer for parameters of model
        loss_fn: a function that takes batch_output and batch_labels and computes the loss for the batch
        dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches training data
        metrics: (dict) a dictionary of functions that compute a metric using the output and labels of each batch
        params: (Params) hyperparameters
        num_steps: (int) number of batches to train on, each of size params.batch_size
    """

    # set model to training mode
    model.train()

    # summary for current training loop and a running average object for loss
    summ = []
    loss_avg = utils.RunningAverage()

    # Use tqdm for progress bar
    # with tqdm(total=len(dataloader)) as t:
    with tqdm(total=params.dict['num_batches_per_epoch']) as t:
        for i, (features, survival) in zip(range(params.dict['num_batches_per_epoch']), dataloader):
            # labels_batch = survival[:, 1]
            train_batch, labels_batch = torch.from_numpy(features).float(), torch.from_numpy(survival[:, 1]).float()
            # move to GPU if available
            if params.cuda:
                train_batch, labels_batch = train_batch.cuda(async=True), labels_batch.cuda(async=True)
            # convert to torch Variables
            # train_batch, labels_batch = Variable(train_batch), Variable(labels_batch)

            # compute model output and loss
            output_batch = model(train_batch)
            loss = loss_fn(labels_batch, output_batch)
            if(torch.isnan(loss)):
                import ipdb
                ipdb.set_trace()
            # import ipdb
            # ipdb.set_trace()

            # clear previous gradients, compute gradients of all variables wrt loss
            optimizer.zero_grad()
            loss.backward()

            # performs updates using calculated gradients
            optimizer.step()

            # output_batch1 = model(train_batch)
            # if(torch.any(torch.isnan(output_batch1))):
            #     import ipdb
            #     ipdb.set_trace()

            # Evaluate summaries only once in a while
            if i % params.save_summary_steps == 0:
                # extract data from torch Variable, move to cpu, convert to numpy arrays
                output_batch = output_batch.data.detach().cpu().numpy()
                labels_batch = labels_batch.data.detach().cpu().numpy()
                # c_index = metrics.concordance_metric(output_batch, survival)

                # compute all metrics on this batch
                # summary_batch = {'c_index': c_index}
                summary_batch = {metric: metrics[metric](output_batch, survival)
                                 for metric in metrics}
                summary_batch['loss'] = loss.item()
                summ.append(summary_batch)

            # update the average loss
            loss_avg.update(loss.item())

            t.set_postfix(loss='{:05.3f}'.format(loss_avg()))
            t.update()

    # compute mean of all metrics in summary
    metrics_mean = {metric: np.mean([x[metric] for x in summ]) for metric in summ[0]}
    metrics_string = " ; ".join("{}: {:05.3f}".format(k, v) for k, v in metrics_mean.items())
    logging.info("- Train metrics: " + metrics_string)


def train_and_evaluate(model, train_dataloader, val_dataloader, optimizer, loss_fn, metrics, params, model_dir,
                       restore_file=None):
    """Train the model and evaluate every epoch.
    Args:
        model: (torch.nn.Module) the neural network
        train_dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches training data
        val_dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches validation data
        optimizer: (torch.optim) optimizer for parameters of model
        loss_fn: a function that takes batch_output and batch_labels and computes the loss for the batch
        metrics: (dict) a dictionary of functions that compute a metric using the output and labels of each batch
        params: (Params) hyperparameters
        model_dir: (string) directory containing config, weights and log
        restore_file: (string) optional- name of file to restore from (without its extension .pth.tar)
    """
    # reload weights from restore_file if specified
    if restore_file is not None:
        restore_path = os.path.join(args.model_dir, args.restore_file + '.pth.tar')
        logging.info("Restoring parameters from {}".format(restore_path))
        utils.load_checkpoint(restore_path, model, optimizer)

    best_val_acc = 0.5  # for cindex

    for epoch in range(params.num_epochs):
        # Run one epoch
        logging.info("Epoch {}/{}".format(epoch + 1, params.num_epochs))

        # compute number of batches in one epoch (one full pass over the training set)
        train(model, optimizer, loss_fn, train_dataloader, metrics, params)

        # Evaluate for one epoch on validation set
        # val_metrics = evaluate(model, loss_fn, val_dataloader, metrics, params)
        val_metrics = {'c_index': 0.5}

        val_acc = val_metrics['c_index']
        is_best = val_acc > best_val_acc

        # Save weights
        utils.save_checkpoint({'epoch': epoch + 1,
                               'state_dict': model.state_dict(),
                               'optim_dict': optimizer.state_dict()},
                              is_best=is_best,
                              checkpoint=model_dir)

        # If best_eval, best_save_path
        if is_best:
            logging.info("- Found new best cindex")
            best_val_acc = val_acc

            # Save best val metrics in a json file in the model directory
            best_json_path = os.path.join(model_dir, "metrics_val_best_weights.json")
            utils.save_dict_to_json(val_metrics, best_json_path)

        # Save latest val metrics in a json file in the model directory
        last_json_path = os.path.join(model_dir, "metrics_val_last_weights.json")
        utils.save_dict_to_json(val_metrics, last_json_path)
        # return the trained model


if __name__ == '__main__':

    # Load the parameters from json file
    args = parser.parse_args()
    json_path = os.path.join(args.model_dir, 'params.json')
    assert os.path.isfile(json_path), "No json configuration file found at {}".format(json_path)
    params = utils.Params(json_path)
    # print(type(params))
    # print(params)

    # use GPU if available
    params.cuda = torch.cuda.is_available()

    # Set the random seed for reproducible experiments
    torch.manual_seed(230)
    # if params.cuda:
    # torch.cuda.manual_seed(230)

    # Set the logger
    utils.set_logger(os.path.join(args.model_dir, 'train.log'))

    # Create the input data pipeline
    logging.info("Loading the datasets...")

    # fetch dataloaders
    dataloaders = data_generator.fetch_dataloader(['train', 'val'], args.data_dir, params)
    train_steps_gen, train_input_size, train_dl = dataloaders['train']
    _, _, val_dl = dataloaders['val']
    input_size = train_input_size
    params.dict['num_batches_per_epoch'] = train_steps_gen
    logging.info("- done.")

    # Define the model and optimizer
    # model = net.FCN(params).cuda() if params.cuda else net.FCN(params)
    model = net.NeuralNet(input_size, params.hidden_size, 1).cuda() if params.cuda else net.NeuralNet(input_size, params.hidden_size, 1)
    optimizer = optim.Adam(model.parameters(), lr=params.learning_rate, weight_decay=1e-3)

    # fetch loss function and metrics
    loss_fn = net.negative_log_partial_likelihood
    metrics = net.metrics

    # Train the model
    logging.info("Starting training for {} epoch(s)".format(params.num_epochs))
    train_and_evaluate(model, train_dl, val_dl, optimizer, loss_fn, metrics, params, args.model_dir,
                       args.restore_file)
