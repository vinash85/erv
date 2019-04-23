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
parser.add_argument('--hyper_param', default='',
                    help="support string for setting parameter from command line e.g.\"params.input_indices=range(50)\"")
parser.add_argument('--type_file', default="val",
                    help="String for files from dataset for validation list will be executed")


def evaluate_attention(models, dataloader, metrics, params, validation_file=None, writer=None, epoch=None, index=None, tsne=0):
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
    for inx in range(len(models)):
        models[inx].eval()

    num_batches_per_epoch, _, _, dataloader = dataloader

    # summary for current eval loop
    summ = []

    # compute metrics over the dataset
    # print(num_batches_per_epoch)
    # import ipdb
    # ipdb.set_trace()

    for i, (features, all_labels, curr_response) in zip(range(num_batches_per_epoch), dataloader):
        # labels_san_survival = np.take(all_labels, params.survival_indices + params.continuous_phenotype_indices + params.binary_phentoype_indices, axis=1).astype(float)
        labels_san_survival = all_labels
        data_batch, labels_batch = torch.from_numpy(features).float(), torch.from_numpy(labels_san_survival).float()
        if(i == 0):
            response = curr_response
        else:
            response = np.concatenate([response, curr_response], 0)
        # move to GPU if available
        if params.cuda:
            data_batch, labels_batch = data_batch.cuda(non_blocking=True), labels_batch.cuda(non_blocking=True)
        # fetch the next evaluation batch
        embedding_input = data_batch[:, params.embedding_indices]
        attention_input = data_batch[:, params.attention_indices]

        embedding_batch = models[0](embedding_input)
        attention_mat = models[1](attention_input)
        transformed_batch = net.feature_attention(attention_mat, embedding_batch)
        output_batch = models[2](transformed_batch)

        # compute model output
        loss = net.update_loss_parameters_vectorized(labels_batch, output_batch, models, None, params, [0, 0, 0])

        # extract data from torch Variable, move to cpu, convert to numpy arrays
        output_batch = output_batch.data.cpu().numpy()
        embedding_batch = embedding_batch.data.cpu().numpy()
        transformed_batch = transformed_batch.data.cpu().numpy()

        # if validation_file and writer is not None:
        output_and_predictions = np.concatenate([labels_san_survival, output_batch, embedding_batch, transformed_batch], 1)
        if i == 0:
            predictions = output_and_predictions
        else:
            predictions = np.concatenate([predictions, output_and_predictions], 0)
        # labels_batch = labels_batch.data.cpu().numpy()

        # compute all metrics on this batch
        summary_batch = {dd[0] + "_" + dd[1]: metrics[dd[1]](output_batch[:, dd[2]], labels_san_survival[:, dd[3]: (dd[3] + 2)])
                         if dd[1] == 'c_index' else
                         metrics[dd[1]](output_batch[:, dd[2]], labels_san_survival[:, dd[3]])
                         for inx, dd in enumerate(params.metrics)}  # TODO ugly solution, when more metrics change it!!

        summary_batch['loss'] = loss
        summary_batch['negative_loss'] = -loss
        summ.append(summary_batch)

    # compute mean of all metrics in summary
    # print(summ)
    metrics_mean = {metric: net.mean_na([x[metric] for x in summ]) for metric in summ[0]}
    metrics_string = " ; ".join("{}: {:05.3f}".format(k, v) for k, v in metrics_mean.items())
    logging.info("- Eval metrics : " + metrics_string)
    if validation_file:
        predictions_df = pd.DataFrame(predictions)
        header_outputs = [params.header[ii] + ".output" for ii in list(range(0, len(params.survival_indices), 2)) + list(range(len(params.survival_indices), len(params.header)))]
        embedding_header = ["embedding" + str(inx) for inx in range(params.embedding_size)]
        transformed_header = ["transformed" + str(inx) for inx in range(transformed_batch.shape[1])]
        all_header = params.header.tolist() + header_outputs + embedding_header + transformed_header
        predictions_df.columns = all_header
        predictions_df.to_csv(validation_file, sep='\t', index=False)

    if tsne:
          # will work with epoch
        if epoch % params.embedding_log == 0:
            split_indices = np.cumsum([labels_san_survival.shape[1], output_batch.shape[1],
                                       embedding_batch.shape[1]])
            _, output_mat, embedding_mat, transformed_mat = \
                np.split(predictions, split_indices, axis=1)

            response = [tuple(row) for row in response]
            label_img = np.apply_along_axis(utils.reshape_toimage, 1, output_mat)
            writer.add_embedding(
                embedding_mat,
                metadata=response,
                label_img=torch.from_numpy(label_img).float(),
                global_step=epoch,
                metadata_header=params.metadata_header,
                tag='val_' + 'embedding' + str(index))

            writer.add_embedding(
                transformed_mat,
                metadata=response,
                label_img=torch.from_numpy(label_img).float(),
                global_step=epoch,
                metadata_header=params.metadata_header,
                tag='val_' + 'transformed' + str(index))

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
    exec(args.hyper_param)

    # Set the random seed for reproducible experiments
    torch.manual_seed(230)
    if params.cuda:
        torch.cuda.manual_seed(230)

    # Get the logger
    utils.set_logger(os.path.join(args.model_dir, 'evaluate.log'))

    # Create the input data pipeline
    # params.survival_indices = eval(params.survival_indices)
    # params.continuous_phenotype_indices = eval(params.continuous_phenotype_indices)
    # params.binary_phentoype_indices = eval(params.binary_phentoype_indices)

    # params.loss_excluded_from_training = eval(params.loss_excluded_from_training)
    # params.metrics = eval(params.metrics)

    logging.info("Creating the dataset...")
    # net.tracer()
    # fetch dataloaders
    type_file = args.type_file
    # fetch dataloaders
    datasets = data_generator.fetch_dataloader_list(args.prefix,
                                                    [type_file], args.data_dir, params, shuffle=False)

    params = net.create_lossfns_mask(params)
    _, _, params.header, _ = datasets[0][0][type_file]
    params.input_size = len(params.embedding_indices)
    params.attention_input_size = len(params.attention_indices)
    params = net.define_metrics(params)
    logging.info("- done.")

    # Define the model and optimizer
    modelClasses = [net.EmbeddingNet, net.AttentionEncoder, net.outputLayer]
    models = [modelClass(params).cuda() if params.cuda else modelClass(params) for modelClass in modelClasses]

    # Define the model

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
    utils.load_checkpoint_attn(restore_path, models)

    metrics = net.metrics

    # Evaluate for one epoch on validation set

    data_dirs = pd.read_csv(args.data_dir, sep="\t")
    data_dirs = [row['data_dir'] for index, row in data_dirs.iterrows()]
    val_metrics_all = [evaluate_attention(models, dataset[0][type_file], metrics, params, validation_file=data_dir + "/" + type_file + "_prediction.csv") for data_dir, dataset in zip(data_dirs, datasets)]
    val_metrics = {metric: eval(params.aggregate)([x[metric] for x in val_metrics_all]) for metric in val_metrics_all[0]}
    save_path = os.path.join("{}.json".format(restore_path))
    utils.save_dict_to_json(val_metrics, save_path)
