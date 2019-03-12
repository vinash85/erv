"""Evaluates the model"""

import argparse
import logging
import os

import numpy as np
import torch
from torch.autograd import Variable
import utils
import pandas as pd
import model.net as net
import model.data_generator as data_generator


parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', default='data/64x64_SIGNS', help="Directory containing the dataset")
parser.add_argument('--model_dir', default='experiments/base_model', help="Directory containing params.json")
parser.add_argument('--restore_file', default='best', help="name of the file in --model_dir \
                     containing weights to load")
parser.add_argument('--prefix', default='',
                    help="Prefix of dataset files  \n \
                    (e.g. prefix=\"tcga\" implies input files are \n \
                    tcga_ssgsea_[train,test,val].txt, \n \
                    tcga_phenotype_[train,test,val].txt )")


def evaluate(embedding_model, outputs, dataloader, metrics, params, validation_file=None):
    """Evaluate the model on `num_steps` batches.

    Args:
        model: (torch.nn.Module) the neural network
        dataloader: (DataLoader) a torch.utils.data.DataLoader object that fetches data
        metrics: (dict) a dictionary of functions that compute a metric using the output and labels of each batch
        params: (Params) hyperparameters
        num_steps: (int) number of batches to train on, each of size params.batch_size
    """

    # set model to evaluation mode
    # print(embedding_model.state_dict())
    embedding_model.eval()
    outputs.eval()
    predictions = np.array([])

    num_batches_per_epoch, _, dataloader = dataloader

    # summary for current eval loop
    summ = []

    # compute metrics over the dataset
    # print(num_batches_per_epoch)
    # import ipdb
    # ipdb.set_trace()

    for i, (features, all_labels) in zip(range(num_batches_per_epoch), dataloader):
        labels_san_survival = np.take(all_labels, params.survival_indices + params.continuous_phenotype_indices + params.binary_phentoype_indices, axis=1)
        data_batch, labels_batch = torch.from_numpy(features).float(), torch.from_numpy(labels_san_survival).float()
        # move to GPU if available
        if params.cuda:
            data_batch, labels_batch = data_batch.cuda(non_blocking=True), labels_batch.cuda(non_blocking=True)
        # fetch the next evaluation batch

        # compute model output
        embedding_batch = embedding_model(data_batch)
        output_batch = outputs(embedding_batch)
        loss = net.update_loss_parameters(
            labels_batch, output_batch, embedding_model, outputs, None, None, params, [0, 0])
        # extract data from torch Variable, move to cpu, convert to numpy arrays
        output_batch = output_batch.data.cpu().numpy()

        if validation_file:
            output_and_predictions = np.concatenate([labels_san_survival, output_batch], 1)
            if i == 0:
                predictions = output_and_predictions
            else:
                predictions = np.concatenate([predictions, output_and_predictions], 0)
        # labels_batch = labels_batch.data.cpu().numpy()

        # compute all metrics on this batch
        summary_batch = {dd[0]: metrics[dd[0]](output_batch[:, dd[1]], labels_san_survival[:, dd[2]: (dd[2] + 2)])
                         if dd[0] == 'c_index' else
                         metrics[dd[0]](output_batch[:, dd[1]], labels_san_survival[:, dd[2]])
                         for dd in params.metrics}

        summary_batch['loss'] = loss
        summary_batch['negative_loss'] = -loss
        summ.append(summary_batch)

    # compute mean of all metrics in summary
    # print(summ)
    metrics_mean = {metric: net.mean_na([x[metric] for x in summ]) for metric in summ[0]}
    metrics_string = " ; ".join("{}: {:05.3f}".format(k, v) for k, v in metrics_mean.items())
    logging.info("- Eval metrics : " + metrics_string)
    if validation_file:
        predictions = pd.DataFrame(predictions)
        predictions.to_csv(validation_file, sep='\t', index=False)

    return metrics_mean


if __name__ == '__main__':
    """
        Evaluate the model on the test set.
    """
    # Load the parameters
    args = parser.parse_args()
    json_path = os.path.join(args.model_dir, 'params.json')
    assert os.path.isfile(json_path), "No json configuration file found at {}".format(json_path)
    params = utils.Params(json_path)

    # use GPU if available
    params.cuda = torch.cuda.is_available()     # use GPU is available

    # Set the random seed for reproducible experiments
    torch.manual_seed(230)
    if params.cuda:
        torch.cuda.manual_seed(230)

    # Get the logger
    utils.set_logger(os.path.join(args.model_dir, 'evaluate.log'))

    # Create the input data pipeline
    logging.info("Creating the dataset...")

    # fetch dataloaders
    type_file = 'val'  # or 'test'

    # fetch dataloaders
    datasets = data_generator.fetch_dataloader_list(args.prefix,
                                                    [type_file], args.data_dir, params, shuffle=False)
    _, input_size, _ = datasets[0][0][type_file]

    logging.info("- done.")
    params.survival_indices = np.asarray(eval(params.survival_indices), dtype=np.int)
    params.continuous_phenotype_indices = np.asarray(eval(params.continuous_phenotype_indices), dtype=np.int)
    params.binary_phentoype_indices = np.asarray(eval(params.binary_phentoype_indices), dtype=np.int)

    params.loss_excluded_from_training = np.asarray(eval(params.loss_excluded_from_training), dtype=np.int)

    params.loss_fns, params.mask, linear_output_size, binary_output_size = net.create_lossfns_mask(params)

    # Define the model
    embedding_model = net.EmbeddingNet(
        net.ConvolutionBlock, input_size, out_channels_list=params.out_channels_list, FC_size_list=params.FC_size_list, embedding_size=params.embedding_size, kernel_sizes=params.kernel_sizes, strides=params.strides, dropout_rate=params.dropout_rate)

    outputs = net.outputLayer_simple(params.embedding_size, linear_output_size=linear_output_size,
                                     binary_output_size=binary_output_size)

    if params.cuda:
        # model = model.cuda()
        embedding_model = embedding_model.cuda()
        outputs = outputs.cuda()

    metrics = net.metrics

    logging.info("Starting evaluation")

    # Reload weights from the saved file
    # print(os.path.join(args.model_dir, args.restore_file + '.pth.tar'))
    if os.path.isfile(args.restore_file):
        restore_path = args.restore_file
    else:
        restore_path = os.path.join(
            args.model_dir, args.restore_file + '.pth.tar')
    logging.info("Restoring parameters from {}".format(restore_path))
    utils.load_checkpoint(restore_path, embedding_model, outputs)

    metrics = net.metrics

    # Evaluate for one epoch on validation set

    data_dirs = pd.read_csv(args.data_dir, sep="\t")
    data_dirs = [row['data_dir'] for index, row in data_dirs.iterrows()]
    val_metrics_all = [evaluate(embedding_model, outputs, dataloader[0][type_file], metrics, params, validation_file=data_dir + "/" + type_file + "_prediction.csv") for data_dir, dataloader in zip(data_dirs, datasets)]
    val_metrics = {metric: eval(params.aggregate)([x[metric] for x in val_metrics_all]) for metric in val_metrics_all[0]}
    save_path = os.path.join("{}.json".format(restore_path))
    utils.save_dict_to_json(val_metrics, save_path)
