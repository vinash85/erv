"""Defines the neural network, loss function and metrics"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from lifelines.utils import concordance_index
import math
from sklearn import metrics as smetrics
import r2python

# 3x3 convolution


def conv1d(in_channels, out_channels, kernel_size=5, stride=2):
    padding_size = int((kernel_size - 1) / 2)
    # print(padding_size)
    # print(in_channels)
    return nn.Conv1d(in_channels, out_channels, kernel_size=kernel_size,
                     stride=stride, padding=padding_size, bias=False)

# Residual block


class ConvolutionBlock(nn.Module):

    def __init__(self, in_channels, out_channels, kernel_size=5, stride=2):
        super(ConvolutionBlock, self).__init__()
        # print(in_channels)
        # print(out_channels)
        self.conv1 = conv1d(in_channels, out_channels, kernel_size, stride=1)
        self.bn1 = nn.BatchNorm1d(out_channels)
        # self.relu = nn.Residual(inplace=True)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool1d(kernel_size=kernel_size, stride=stride)

    def forward(self, x):
        residual = x
        out = self.conv1(residual)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.maxpool(out)
        return out


# tempNet
class tempNet(nn.Module):

    def __init__(self, block, input_size, out_channels_list=[32, 32, 32],
                 embedding_size=32, blocks=1, kernel_sizes=[5], strides=[2]):
        super(tempNet, self).__init__()
        self.in_channels = 1
        self.block = block
        if len(kernel_sizes) == 1:
            self.kernel_sizes = [kernel_sizes] * len(out_channels_list)
        if len(strides) == 1:
            self.strides = [strides] * len(out_channels_list)

        self.lay1 = self.make_layer(block, out_channels_list[
                                    0], kernel_sizes=self.kernel_sizes, strides=self.strides[0])
        # self.embedding_size = embedding_size
        # self.blocks = blocks
        # self.kernel_size = kernel_size
        # self.strides = strides
       # (1, 64, 5, stride=2

    def make_layer(self, block, out_channels, kernel_sizes, strides):
        layers = []
        layers.append(block(self.in_channels, out_channels,
                            kernel_sizes[0], stride))
        self.in_channels = out_channels
        for i in range(1, blocks):
            layers.append(
                block(out_channels, out_channels, kernel_size, stride=1))
        return nn.Sequential(*layers)

    def forward(self, x):
        out = x
        # out = self.make_layer(self.block,  self.out_channels_list[0], kernel_size = self.kernel_size, stride = self.strides[0], blocks=self.blocks)(out)
        out = self.lay1(out)
        print("here")
        print(out)
        return out


# EmbeddingNet


# EmbeddingNet
class EmbeddingNet(nn.Module):

    def __init__(self, block, input_size, out_channels_list,
                 embedding_size=32, kernel_sizes=[5], strides=[2], FC_size_list=[32],
                 dropout_rate=0.1):
        super(EmbeddingNet, self).__init__()
        self.in_channels = 1
        self.kernel_sizes = kernel_sizes
        self.strides = strides
        self.dropout_rate = dropout_rate
        if len(kernel_sizes) == 1:
            self.kernel_sizes = [kernel_sizes[0]] * len(out_channels_list)
        if len(strides) == 1:
            self.strides = [strides[0]] * len(out_channels_list)
        self.output_size = input_size
        print("initial output size")
        print(self.output_size)
        self.layers_block1 = self.make_layers(
            block, out_channels_list, kernel_sizes=self.kernel_sizes, strides=self.strides)
        # output_size is updated
        print("final convolution layer output size")
        self.fc_output_size = self.output_size * out_channels_list[-1]

        # self.fc2 = nn.Linear(self.fc_input_size, 2 * embedding_size)
        # self.fc2 = nn.Linear(self.fc_input_size, 2 * embedding_size)
        # print(self.fc_input_size)
        # self.fc3 = nn.Linear(2 * embedding_size, embedding_size)

        print("initial fully connected size")
        print(self.fc_output_size)
        self.FC_block1 = self.make_layers_FC(
            FullConnectedBlock, FC_size_list, dropout_rate)
        print("final output size")
        self.fc3 = nn.Linear(self.fc_output_size, embedding_size)

    def make_layers_FC(self, block, FC_size_list, dropout_rate):
        layers = []
        num_layers = len(FC_size_list)
        for i in range(0, num_layers):
            layers.append(block(self.fc_output_size, FC_size_list[
                          i], dropout_rate))
            self.fc_output_size = FC_size_list[i]
            print(self.fc_output_size)
        return nn.Sequential(*layers)

    def make_layers(self, block, out_channels_list, kernel_sizes, strides):
        layers = []
        num_layers = len(out_channels_list)
        for i in range(0, num_layers):
            layers.append(block(self.in_channels, out_channels_list[
                          i], kernel_sizes[i], stride=strides[i]))
            self.in_channels = out_channels_list[i]
            padding_size = (self.kernel_sizes[i] - 1) / 2
            convolution_stride = 1

            self.output_size = int(((self.output_size + 2 * padding_size - kernel_sizes[
                                   i]) / convolution_stride) + 1)  # convolution layer output_size
            # maxpool output_size
            self.output_size = int((
                (self.output_size - kernel_sizes[i]) / strides[i]) + 1)
            print(self.output_size)
        return nn.Sequential(*layers)

    def forward(self, x):
        # reshape numbatch * num_dim to numbatch * num_in_channel * num_dim
        out = x.view(x.size(0), 1, -1)
        out = self.layers_block1(out)
        temp = out.size()
        out = out.view(out.size(0), -1)
        out = self.FC_block1(out)
        out = self.fc3(out)
        return out


class outputLayer(nn.Module):

    def __init__(self, embedding_size, linear_output_size=32, binary_output_size=32):
        super(outputLayer, self).__init__()
        self.linear_output_size = linear_output_size
        self.binary_output_size = binary_output_size
        if linear_output_size > 0:
            self.linear1 = nn.Linear(embedding_size, linear_output_size)
            self.dense1_bn = nn.BatchNorm1d(1)
        if binary_output_size > 0:
            self.linear2 = nn.Linear(embedding_size, binary_output_size)
        self.sigmoid = nn.Sigmoid()
        # self.softmax = nn.LogSoftmax()

    def forward(self, x):
        if(self.linear_output_size > 0):
            linear1_out = self.linear1(x)
            # print(linear1_out.shape)
            # add a batch normalization for survival prediction
            survival_out = self.dense1_bn(linear1_out[:, 0:1])

        if(self.binary_output_size > 0):
            linear2_out = self.linear2(x)
            # binary_output = self.sigmoid(linear2_out)
            binary_output = linear2_out  # sigmoid not required with BCEWithLogitsLoss
            binary_output = self.sigmoid(binary_output)

        if(self.linear_output_size == 0):
            out = binary_output
        elif(self.binary_output_size == 0):
            out = torch.cat((survival_out, linear1_out[:, 1:]), 1)
        else:
            out = torch.cat((survival_out, linear1_out[:, 1:], binary_output), 1)
        # print(out.shape)
        return out

# FC block


class FullConnectedBlock(nn.Module):

    def __init__(self, in_channels, out_channels, dropout_rate):
        super(FullConnectedBlock, self).__init__()
        self.fc = nn.Linear(in_channels, out_channels)
        self.bn1 = nn.BatchNorm1d(out_channels)
        self.dropout_rate = dropout_rate
        self.relu = nn.ReLU(inplace=True)

    def forward(self, x):
        residual = x
        out = self.fc(residual)
        out = self.bn1(out)
        out = self.relu(out)
        out = F.dropout(out, p=self.dropout_rate, training=self.training)
        return out


# EmbeddingNet
class EmbeddingNet_FC(nn.Module):

    def __init__(self, block, input_size, out_channels_list,
                 embedding_size=32,
                 dropout_rate=0.1):
        super(EmbeddingNet_FC, self).__init__()
        self.in_channels = input_size
        self.output_size = input_size
        self.dropout_rate = dropout_rate
        print("initial output size")
        print(self.output_size)
        self.layers_block1 = self.make_layers_FC(
            block, out_channels_list, dropout_rate)
        # output_size is updated
        print("final output size")
        print(self.output_size)
        self.fc3 = nn.Linear(self.in_channels, embedding_size)

    def make_layers_FC(self, block, out_channels_list, dropout_rate):
        layers = []
        num_layers = len(out_channels_list)
        for i in range(0, num_layers):
            layers.append(block(self.in_channels, out_channels_list[
                          i], dropout_rate))
            self.in_channels = out_channels_list[i]
        return nn.Sequential(*layers)

    def forward(self, x):
        # reshape numbatch * num_dim to numbatch * num_in_channel * num_dim
        # out = x.view(x.size(0), 1, -1)
        out = x
        out = self.layers_block1(out)
        out = self.fc3(out)
        return out


class embedding_conv(nn.Module):

    def __init__(self, input_size, num_classes=1):
        super(ConvNet1D, self).__init__()
        layer1_size = math.floor(input_size / 2)
        output_size = math.floor(layer1_size / 2)
        self.layer1 = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=5, stride=1, padding=2),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2))
        self.layer2 = nn.Sequential(
            nn.Conv1d(32, 32, kernel_size=5, stride=1, padding=2),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2))
        # self.output_size = math_ceiling
        self.fc = nn.Linear(output_size * 32, num_classes)
        self.dense1_bn = nn.BatchNorm1d(num_classes)


class oneLayerNet(nn.Module):

    def __init__(self, input_size, num_classes):
        super(oneLayerNet, self).__init__()

    def forward(self, x):
        out = nn.Linear(self.input_size, self.num_classes)(x)
        out = nn.ReLU()(out)
        out = nn.BatchNorm1d(self.num_classes)(out)
        return out


class ImmuneEmbed(nn.Module):

    def __init__(self, input_size, hidden_layers):
        super(ImmuneEmbed, self).__init__()

    def forward(self, x):
        out = nn.Linear(self.input_size, self.hidden_layer[0])(x)
        out = nn.ReLU()(out)
        out = nn.BatchNorm1d(self.hidden_layer[0])(out)
        for layer_size in self.hidden_layers:
            out = nn.Linear(layer_size)(out)
            out = nn.ReLU()(out)
            out = nn.BatchNorm1d(layer_size)(out)
        return out


# Fully connected neural network with one hidden layer

class NeuralNet(nn.Module):

    def __init__(self, input_size, hidden_size, num_classes=1):
        super(NeuralNet, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        # batchnorm takes both channels or hidden size as input
        self.dense1_bn = nn.BatchNorm1d(hidden_size)
        self.fc2 = nn.Linear(hidden_size, num_classes)
        self.dense2_bn = nn.BatchNorm1d(num_classes)
        # self.dense2_bn = nn.BatchNorm1d(num_classes)

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.dense1_bn(out)
        out = self.fc2(out)
        # print(out.shape)
        out = self.dense2_bn(out)
        return out

# 1D Convolutional neural network (two convolutional layers)


class ConvNet1D(nn.Module):

    def __init__(self, input_size, num_classes=1):
        super(ConvNet1D, self).__init__()
        layer1_size = math.floor(input_size / 2)
        output_size = math.floor(layer1_size / 2)
        self.layer1 = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=5, stride=1, padding=2),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2))
        self.layer2 = nn.Sequential(
            nn.Conv1d(32, 32, kernel_size=5, stride=1, padding=2),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2))
        # self.output_size = math_ceiling
        self.fc = nn.Linear(output_size * 32, num_classes)
        self.dense1_bn = nn.BatchNorm1d(num_classes)

    def forward(self, x):
        gene_len, batch_size = x.size()  # x is gene_len x batch_size
        # print(x.shape)
        # Turn (gene_len x batch_size) into (batch_size x input_size x gene_len) for CNN, where input_size = 1
        # x = x.transpose(0, 1)
        # print(x.shape)
        x = x.reshape(batch_size, 1, gene_len)
        # print(x.shape)
        out = self.layer1(x)
        # print(out.shape)
        out = self.layer2(out)
        out = out.reshape(out.size(0), -1)
        out = self.fc(out)
        out = self.dense1_bn(out)
        return out

# Convolutional neural network (two convolutional layers)


class ConvNet(nn.Module):

    def __init__(self, num_classes=10):
        super(ConvNet, self).__init__()
        self.layer1 = nn.Sequential(
            nn.Conv2d(1, 16, kernel_size=5, stride=1, padding=2),
            nn.BatchNorm2d(16),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2))
        self.layer2 = nn.Sequential(
            nn.Conv2d(16, 32, kernel_size=5, stride=1, padding=2),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2))
        self.fc = nn.Linear(7 * 7 * 32, num_classes)

    def forward(self, x):
        out = self.layer1(x)
        out = self.layer2(out)
        out = out.reshape(out.size(0), -1)
        out = self.fc(out)
        return out


def loss_fn(outputs, labels):
    """
    Compute the cross entropy loss given outputs and labels.

    Args:
        outputs: (Variable) dimension batch_size x 6 - output of the model
        labels: (Variable) dimension batch_size, where each element is a value in [0, 1, 2, 3, 4, 5]

    Returns:
        loss (Variable): cross entropy loss for all images in the batch

    Note: you may use a standard loss function from http://pytorch.org/docs/master/nn.html#loss-functions. This example
          demonstrates how you can easily define a custom loss function.
    """
    num_examples = outputs.size()[0]
    return -torch.sum(outputs[range(num_examples), labels]) / num_examples


def checkNaNInf(x):
    torch.isnan(x).any() or x.eq(
        float('inf')).any() or x.eq(float('-inf')).any()


def isfinite(x):
    """
    Quick pytorch test that there are no nan's or infs.

    note: torch now has torch.isnan
    url: https://gist.github.com/wassname/df8bc03e60f81ff081e1895aabe1f519
    """
    not_inf = ((x + 1) != x)
    not_nan = (x == x)
    return not_inf & not_nan


def negative_log_partial_likelihood(survival, risk, debug=False):
    """Return the negative log-partial likelihood of the prediction
    y_true contains the survival time
    risk is the risk output from the neural network
    censor is the vector of inputs that are censored
    regularization is the regularization constant (not used currently in model)

    Uses the torch backend to perform calculations

    Sorts the surv_time by sorted reverse time
    """

    # calculate negative log likelihood from estimated risk\

    # temp = survival
    # temp[:, 1] = risk

    # print(torch.stack([survival[:, 0], risk]))
    _, idx = torch.sort(survival[:, 0], descending=True)
    censor = survival[idx, 1]
    risk = risk[idx]
    # print(temp[idx, :])
    # print(survival[idx, :])
    epsilon = 0.00001
    max_value = 10
    alpha = 0.1
    risk = torch.reshape(risk, [-1])  # flatten
    # risk = risk - torch.mean(risk)
    # risk[risk > max_value] = max_value
    hazard_ratio = torch.exp(risk)

    # cumsum on sorted surv time accounts for concordance
    log_risk = torch.log(torch.cumsum(hazard_ratio, dim=0) + epsilon)
    log_risk = torch.reshape(log_risk, [-1])
    uncensored_likelihood = risk - log_risk

    # apply censor mask: 1 - dead, 0 - censor
    censored_likelihood = uncensored_likelihood * censor
    num_observed_events = torch.sum(censor)
    neg_likelihood = - torch.sum(censored_likelihood) / \
        num_observed_events
    # if(not isfinite(neg_likelihood) or debug):
    #     import ipdb
    #     ipdb.set_trace()

    # print(type(neg_likelihood))

    # if (neg_likelihood != neg_likelihood):
    #     print(neg_likelihood)
    #     print(censor[np.isnan(censor)])
    #     print(risk[np.isnan(risk)])
    #     raise ValueError("nan found")

    # return torch.exp(neg_likelihood)
    return neg_likelihood


def accuracy(outputs, labels):
    """
    Compute the accuracy, given the outputs and labels for all images.

    Args:
        outputs: (np.ndarray) dimension batch_size x 6 - log softmax output of the model
        labels: (np.ndarray) dimension batch_size, where each element is a value in [0, 1, 2, 3, 4, 5]

    Returns: (float) accuracy in [0,1]
    """
    outputs = np.argmax(outputs, axis=1)
    return np.sum(outputs == labels) / float(labels.size)


def negative_log_partial_likelihood_loss(risk, censor, debug=False):
    '''
    function to change the order of censor and risk
    '''

    return negative_log_partial_likelihood(censor, risk, debug)


def c_index(predicted_risk, survival):
    # calculate the concordance index
    ci = 0  # just to know that concordance index cannot be estimated
    # print(r2python.cbind(np.reshape(predicted_risk, (-1, 1)), survival))

    na_inx = ~(np.isnan(survival[:, 0]) | np.isnan(survival[:, 1]) | np.isnan(predicted_risk))
    predicted_risk, survival = predicted_risk[na_inx], survival[na_inx]
    if len(predicted_risk) > 0:
        survival_time, censor = survival[:, 0], survival[:, 1]
        epsilon = 0.001
        partial_hazard = np.exp(-(predicted_risk + epsilon))
        censor = censor.astype(int)
        ci = concordance_index(survival_time, partial_hazard, censor)
    return ci


# maintain all metrics required in this dictionary- these are used in the
# training and evaluation loops
def calculate_auc(pred, y):
    # print(y)
    # print("pred")
    # print(pred)
    na_inx = ~(isnan(pred) | isnan(y))
    pred, y = pred[na_inx], y[na_inx]
    try:
        fpr, tpr, thresholds = smetrics.roc_curve(y, pred, pos_label=1)
        auc = smetrics.auc(fpr, tpr)
    except:
        auc = 0.5
    return auc

metrics = {
    'auc': calculate_auc,
    'c_index': c_index
    # could add more metrics such as accuracy for each token type
}


def isnan(x):
    return x != x


def calculate_loss(labels, net_outputs, loss_fns):
    '''
    define loss function :
    '''

    # total_loss = torch.zeros(1)
    total_loss = 0.

    # ## survival output
    # survival = labels[:,0:2]
    # na_inx = ~( np.isnan(survival[:,1]) | np.isnan(survival[:,1]))
    # survival, net_output = survival[na_inx,:], net_outputs[na_inx,0]
    # if(len(label) > 1 ):
    #     loss_curr = loss_fns[0](net_output, label)
    # else:
    #     loss_curr = 0

    # total_loss = total_loss + loss_curr

    # label, net_output, loss_fn = labels[:,i], net_outputs[:, i], loss_fns[i]

    len_fns = len(loss_fns)
    label_inx = 0

    for i in range(len_fns):
        net_output, loss_fn = net_outputs[:, i], loss_fns[i]
        # print(loss_fn)
        # print(loss_fn.__name__ is 'negative_log_partial_likelihood_loss')
        if hasattr(loss_fn, '__name__'):
            if loss_fn.__name__ is 'negative_log_partial_likelihood_loss':
                label = labels[:, label_inx:(label_inx + 2)]
                label_inx = label_inx + 2
                # print(label.shape)
                na_inx = ~(isnan(label[:, 0]) | isnan(label[:, 1]))
                label, net_output = label[na_inx, :], net_output[na_inx]
                # print(label)
        else:
            label = labels[:, label_inx]
            label_inx = label_inx + 1
            na_inx = ~isnan(label)
            label, net_output = label[na_inx], net_output[na_inx]

        if(len(label) > 1):
            # in case of survival label is censor; censored data assumed to
            # sorted based on patient event time.
            loss_curr = loss_fn(net_output, label)
        else:
            # loss_curr = torch.zeros(1)
            loss_curr = 0.
        # print(loss_curr.item())
        # manually tries to value response rate
        # if i >= len_fns - 2:
            # loss_curr = loss_curr * 10

        total_loss = total_loss + loss_curr

    return total_loss


def max_na(xx):
    return (0 if len(xx) == 0 else max(xx))


def create_lossfns_mask(params):
    survival_indices = np.asarray(eval(params.survival_indices), dtype=np.int)
    continuous_phenotype_indices = np.asarray(eval(params.continuous_phenotype_indices), dtype=np.int)
    binary_phentoype_indices = np.asarray(eval(params.binary_phentoype_indices), dtype=np.int)

    survival_output_size = int(len(survival_indices) / 2)
    linear_output_size = survival_output_size + len(continuous_phenotype_indices)
    binary_output_size = len(binary_phentoype_indices)
    loss_fns = [negative_log_partial_likelihood_loss] * survival_output_size + \
        [nn.MSELoss()] * len(continuous_phenotype_indices) + \
        [nn.BCELoss()] * len(binary_phentoype_indices)
    # BCEWithLogitsLoss is other option

    # max_index = max(max_na(survival_indices), max_na(continuous_phenotype_indices), max_na(binary_phentoype_indices))

    survival_indices_new, continuous_indices_new, binary_indices_new = survival_indices, continuous_phenotype_indices, binary_phentoype_indices
    print(survival_indices_new)

    # for sur in survival_indices:
    #     print(sur)
    #     print(survival_indices_new > sur)
    #     survival_indices_new[survival_indices_new > sur] = survival_indices_new[survival_indices_new > sur] - 1
    #     continuous_indices_new[continuous_indices_new > sur] = continuous_indices_new[continuous_indices_new > sur] - 1
    #     binary_indices_new[binary_indices_new > sur] = binary_indices_new[binary_indices_new > sur] - 1
    mask = np.concatenate([survival_indices_new, continuous_indices_new, binary_indices_new])
    # print(survival_indices_new)
    print(mask)

    return loss_fns, mask, linear_output_size, binary_output_size


def regularized_loss(model, params, index="all"):

    loss_curr = 0.

    for name, parameters in model.named_parameters():
        if name.endswith('weight') and name.endswith('linear'):
            if index is "all":
                loss_curr = loss_curr + \
                    torch.norm(parameters.view(-1), 1) * params.l1_regularizer \
                    + torch.norm(parameters.view(-1), 2) * params.l2_regularizer
            else:
                loss_curr = loss_curr + \
                    torch.norm(parameters[index].view(-1), 1) * params.l1_regularizer \
                    + torch.norm(parameters[index].view(-1), 2) * params.l2_regularizer
    return loss_curr


def update_loss_parameters(labels, net_outputs, embedding_model, outputs, embedding_optimizer, outputs_optimizer, params, is_train):
    '''
    define loss function and update parameters
    '''

    def update_parameters(loss):
        # clear previous gradients, compute gradients of all variables wrt loss
        embedding_optimizer.zero_grad()
        outputs_optimizer.zero_grad()
        loss.backward(retain_graph=True)
        # performs updates using calculated gradients
        embedding_optimizer.step()
        outputs_optimizer.step()

    total_loss = 0.
    loss_fns = params.loss_fns
    len_fns = len(loss_fns)

    label_inx = 0
    for i in range(len_fns):
        net_output, loss_fn = net_outputs[:, i], loss_fns[i]
        if hasattr(loss_fn, '__name__'):
            if loss_fn.__name__ is 'negative_log_partial_likelihood_loss':
                label = labels[:, label_inx:(label_inx + 2)]
                label_inx = label_inx + 2
                # print(label.shape)
                na_inx = ~(isnan(label[:, 0]) | isnan(label[:, 1]))
                label, net_output = label[na_inx, :], net_output[na_inx]
        else:
            label = labels[:, label_inx]
            label_inx = label_inx + 1
            na_inx = ~isnan(label)
            label, net_output = label[na_inx], net_output[na_inx]

        if(len(label) > 1):
            # in case of survival label is censor; censored data assumed to
            # sorted based on patient event time.
            loss_curr = loss_fn(net_output, label)
            # import ipdb
            # ipdb.set_trace()
            loss_curr = loss_curr + regularized_loss(outputs, params, i)
        else:
            # loss_curr = torch.zeros(1)
            loss_curr = 0.

        if is_train and params.pipeline_optimization:
            update_parameters(loss_curr)

        total_loss = total_loss + loss_curr

    if is_train and not params.pipeline_optimization:
        update_parameters(total_loss)

    return total_loss
