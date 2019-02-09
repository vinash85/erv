"""Train the model"""


import torch
import torch.nn as nn
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
parser.add_argument('--embedding_size', default=32,
                    help="Size of immune embedding (default:8)")
parser.add_argument('--data_dir', default='data/64x64_SIGNS',
                    help="File containing directory containing datasets")
# parser.add_argument('--data_dir_list', default=None,
# help="File contating list of dataset directories data_dirs")
parser.add_argument('--model_dir', default='experiments/base_model',
                    help="Directory containing params.json")
parser.add_argument('--prefix', default='',
                    help="Prefix of dataset files  \n \
                    (e.g. prefix=\"tcga\" implies input files are \n \
                    tcga_ssgsea_[train,test,val].txt, \n \
                    tcga_phenotype_[train,test,val].txt )")
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


def train(embedding_model, outputs, embedding_optimizer, outputs_optimizer, dataloader, metrics, params):
    """Train the model on `num_steps` batches
    Args:
        model: (torch.nn.Module) the neural network
        optimizer: (torch.optim) optimizer for parameters of model
        dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches training data
        metrics: (dict) a dictionary of functions that compute a metric using the output and labels of each batch
        params: (Params) hyperparameters
        num_steps: (int) number of batches to train on, each of size params.batch_size
    """

    # set model to training mode
    embedding_model.train()
    outputs.train()
    num_batches_per_epoch, _, dataloader = dataloader

    # summary for current training loop and a running average object for loss
    summ = []
    loss_avg = utils.RunningAverage()

    # Use tqdm for progress bar
    # with tqdm(total=len(dataloader)) as t:
    with tqdm(total=num_batches_per_epoch) as t:
        for i, (features, all_labels) in zip(range(num_batches_per_epoch), dataloader):
            survival = all_labels[:, 0:2]
            mask = np.ones(all_labels.shape[1], dtype=bool)
            mask[[0, all_labels.shape[1] - 3]] = False
            # labels_san_survival = all_labels[:, 1:]
            labels_san_survival = all_labels[:, mask]
            # print(all_labels.shape)
            train_batch, labels_batch = torch.from_numpy(
                features).float(), torch.from_numpy(labels_san_survival).float()
            # move to GPU if available
            if params.cuda:
                train_batch, labels_batch = train_batch.cuda(
                    non_blocking=True), labels_batch.cuda(non_blocking=True)

            # convert to torch Variables
            # train_batch, labels_batch = Variable(train_batch), Variable(labels_batch)

            # compute model output and loss
            embedding_batch = embedding_model(train_batch)
            # print("embedding_batch")
            # print(embedding_batch.shape)
            output_batch = outputs(embedding_batch)
            loss = net.calculate_loss(
                labels_batch, output_batch, params.loss_fns)
            # params.lambda1 = 0.5
            if params.params_regularization > 0:
                if params.linear_output_size > 0:
                    all_linear1_params = torch.cat([x.view(-1) for x in outputs.linear1.parameters()])
                    l1_regularization = params.params_regularization * torch.norm(all_linear1_params, 2)
                    loss = loss + l1_regularization

                if params.binary_output_size > 0:
                    all_linear2_params = torch.cat([x.view(-1) for x in outputs.linear2.parameters()])
                    l2_regularization = params.params_regularization * torch.norm(all_linear2_params, 2)
                    loss = loss + l2_regularization

            if(torch.isnan(loss)):
                import ipdb
                ipdb.set_trace()
            # import ipdb
            # ipdb.set_trace()

            # clear previous gradients, compute gradients of all variables wrt loss
            embedding_optimizer.zero_grad()
            outputs_optimizer.zero_grad()

            loss.backward()

            # performs updates using calculated gradients
            embedding_optimizer.step()
            outputs_optimizer.step()

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
                # print(output_batch.shape)
                # print(survival.shape)
                # print("line 133")
                # print(survival[0, :])
                summary_batch = {metric: metrics[metric](output_batch[:, 0], survival) if metric == 'c_index' else metrics[metric](output_batch[:, -1], labels_san_survival[:, -1])
                                 for metric in metrics}  # TODO ugly solution, when more metrics change it!!
                # print("line 134")

                summary_batch['loss'] = loss.item()
                summ.append(summary_batch)

            # update the average loss
            loss_avg.update(loss.item())

            t.set_postfix(loss='{:05.3f}'.format(loss_avg()))
            t.update()

    # compute mean of all metrics in summary
    metrics_mean = {metric: np.mean([x[metric]
                                     for x in summ]) for metric in summ[0]}
    metrics_string = " ; ".join("{}: {:05.3f}".format(k, v)
                                for k, v in metrics_mean.items())
    logging.info("- Train metrics: " + metrics_string)


def train_and_evaluate(embedding_model, outputs, datasets, embedding_optimizer, outputs_optimizer, metrics, params, model_dir,
                       restore_file=None):
    """Train the model and evaluate every epoch.
    Args:
        model: (torch.nn.Module) the neural network
        datasets : list of dataloaders, each containing train_dataloader and val_dataloader
        train_dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches training data
        val_dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches validation data
        optimizer: (torch.optim) optimizer for parameters of model
        metrics: (dict) a dictionary of functions that compute a metric using the output and labels of each batch
        params: (Params) hyperparameters
        model_dir: (string) directory containing config, weights and log
        restore_file: (string) optional- name of file to restore from (without its extension .pth.tar)
    """
    # reload weights from restore_file if specified
    if restore_file is not None:
        restore_path = os.path.join(
            args.model_dir, args.restore_file + '.pth.tar')
        logging.info("Restoring parameters from {}".format(restore_path))
        utils.load_checkpoint(restore_path, model, optimizer)

    best_val_acc = 0.5  # for cindex

    for epoch in range(params.num_epochs):
        # Run one epoch
        logging.info("Epoch {}/{}".format(epoch + 1, params.num_epochs))

        # compute number of batches in one epoch (one full pass over the training set)
        for dataloader in datasets:
            train(embedding_model, outputs, embedding_optimizer,
                  outputs_optimizer, dataloader['train'], metrics, params)

        # Evaluate for one epoch on validation set
        val_metrics_all = [evaluate(embedding_model, outputs, dataloader['val'], metrics, params)for dataloader in datasets]
        # val_metrics = []
        val_metrics = {metric: eval(params.metrics)([x[metric] for x in val_metrics_all]) for metric in val_metrics_all[0]}

        # val_metrics = eval(params.metrics)(val_metrics)
        val_acc = val_metrics['c_index']  # use differnt functions

        is_best = val_acc > best_val_acc

        # Save weights
        utils.save_checkpoint({'epoch': epoch + 1,
                               'embedding_state_dict': embedding_model.state_dict(),
                               'outputs': outputs.state_dict(),
                               'embedding_optim_dict': embedding_optimizer.state_dict(),
                               'outputs_optim_dict': outputs_optimizer.state_dict()
                               },
                              is_best=is_best,
                              checkpoint=model_dir)

        # If best_eval, best_save_path
        if is_best:
            logging.info("- Found new best cindex")
            best_val_acc = val_acc

            # Save best val metrics in a json file in the model directory
            best_json_path = os.path.join(
                model_dir, "metrics_val_best_weights.json")
            utils.save_dict_to_json(val_metrics, best_json_path)

        # Save latest val metrics in a json file in the model directory
        last_json_path = os.path.join(
            model_dir, "metrics_val_last_weights.json")
        utils.save_dict_to_json(val_metrics, last_json_path)
        # return the trained model


if __name__ == '__main__':

    # Load the parameters from json file
    args = parser.parse_args()
    json_path = os.path.join(args.model_dir, 'params.json')
    assert os.path.isfile(
        json_path), "No json configuration file found at {}".format(json_path)
    params = utils.Params(json_path)
    if params.loss_fns == 0:
        params.loss_fns = [net.negative_log_partial_likelihood_loss] + [nn.MSELoss()] * (
            params.linear_output_size - 1) + [nn.BCEWithLogitsLoss()] * (params.binary_output_size)

    print(params.loss_fns)
    print(len(params.loss_fns))

    # use GPU if available
    params.cuda = torch.cuda.is_available()
    print(params.cuda)

    # Set the random seed for reproducible experiments
    torch.manual_seed(230)
    # if params.cuda:
    # torch.cuda.manual_seed(230)

    # Set the logger
    utils.set_logger(os.path.join(args.model_dir, 'train.log'))

    # Create the input data pipeline
    logging.info("Loading the datasets...")

    # fetch dataloaders
    datasets = data_generator.fetch_dataloader_list(args.prefix,
                                                    ['train', 'val'], args.data_dir, params)
    _, train_input_size, _ = datasets[0]['train']
    # _, _, val_dl = dataloaders['val']
    # train_dl = dataloaders['train']
    # val_dl = dataloaders['val']
    input_size = train_input_size
    # params.dict['num_batches_per_epoch'] = train_steps_gen
    logging.info("- done.")

    # Define the model and optimizer
    # model = net.FCN(params).cuda() if params.cuda else net.FCN(params)
    # model = net.NeuralNet(input_size, params.hidden_size, 1)
    # embedding_model = net.EmbeddingNet(
    #     net.ConvolutionBlock, input_size, out_channels_list=[32, 32, 32, 32], embedding_size=params.embedding_size, kernel_sizes=[5, 11, 11, 5], strides=[2, 5, 5, 2])
    embedding_model = net.EmbeddingNet(
        net.ConvolutionBlock, input_size, out_channels_list=params.out_channels_list, embedding_size=params.embedding_size, kernel_sizes=params.kernel_sizes, strides=params.strides, dropout_rate=params.dropout_rate)

    outputs = net.outputLayer(params.embedding_size, linear_output_size=params.linear_output_size,
                              binary_output_size=params.binary_output_size)

    if params.cuda:
        # model = model.cuda()
        embedding_model = embedding_model.cuda()
        outputs = outputs.cuda()

    embedding_optimizer = optim.Adam(
        embedding_model.parameters(), lr=params.learning_rate, weight_decay=1e-3)
    outputs_optimizer = optim.Adam(
        outputs.parameters(), lr=params.learning_rate, weight_decay=1e-3)

    # fetch loss function and metrics
    # loss_fn = net.negative_log_partial_likelihood
    metrics = net.metrics

    # Train the model
    logging.info("Starting training for {} epoch(s)".format(params.num_epochs))
    train_and_evaluate(embedding_model, outputs, datasets, embedding_optimizer, outputs_optimizer, metrics, params, args.model_dir,
                       args.restore_file)
