"""Defines the neural network, loss function and metrics"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from lifelines.utils import concordance_index
import math




# 3x3 convolution
def conv1d(in_channels, out_channels, kernel_size=5, stride=2):
    padding_size = int((kernel_size - 1)/2)
    # print(padding_size)
    # print(in_channels)
    return nn.Conv1d(in_channels, out_channels, kernel_size=kernel_size, 
                     stride=stride, padding=padding_size, bias=False)

# Residual block
class ConvolutionBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size = 5, stride=2):
        super(ConvolutionBlock, self).__init__()
        print(in_channels)
        print(out_channels)
        self.conv1 = conv1d(in_channels, out_channels, kernel_size, stride=1)
        self.bn1 = nn.BatchNorm1d(out_channels)
        # self.relu = nn.Residual(inplace=True)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool1d(kernel_size = kernel_size, stride=stride)
        
    def forward(self, x):
        residual = x
        out = self.conv1(residual)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.maxpool(out)
        return out



# tempNet
class tempNet(nn.Module):
    def __init__(self, block, input_size, out_channels_list = [32, 32, 32], 
     embedding_size = 32, blocks = 1, kernel_sizes = [5], strides = [2]):
        super(tempNet, self).__init__()
        self.in_channels = 1
        self.block = block
        if len(kernel_sizes) ==1:
            self.kernel_sizes = [kernel_sizes] * len(out_channels_list)
        if len(strides) ==1:
            self.strides = [strides] * len(out_channels_list)
        
        self.lay1 = self.make_layer(block,  out_channels_list[0], kernel_sizes = self.kernel_sizes, strides = self.strides[0])
        # self.embedding_size = embedding_size
        # self.blocks = blocks
        # self.kernel_size = kernel_size
        # self.strides = strides
       # (1, 64, 5, stride=2 

    def make_layer(self, block, out_channels, kernel_sizes,  strides):
        layers = []
        layers.append(block(self.in_channels, out_channels, kernel_sizes[0], stride))
        self.in_channels = out_channels
        for i in range(1, blocks):
            layers.append(block(out_channels, out_channels, kernel_size, stride=1))
        return nn.Sequential(*layers)
    
    def forward(self, x):
        out = x
        # out = self.make_layer(self.block,  self.out_channels_list[0], kernel_size = self.kernel_size, stride = self.strides[0], blocks=self.blocks)(out)
        out = self.lay1(out)
        print("here")
        print(out)
        return out


# EmbeddingNet
class EmbeddingNet(nn.Module):
    def __init__(self, block, input_size, out_channels_list, 
     embedding_size = 32, kernel_sizes = [5],  strides = [2]):
        super(EmbeddingNet, self).__init__()
        self.in_channels = 1
        if len(kernel_sizes) ==1:
            self.kernel_sizes = [kernel_sizes[0]] * len(out_channels_list)
        if len(strides) == 1:
            self.strides = [strides[0]] * len(out_channels_list)
        self.output_size = input_size
        print(self.output_size)
        self.layers_block1 = self.make_layers(block, out_channels_list, kernel_sizes = self.kernel_sizes, strides = self.strides)
        ## output_size is updated 
        print(self.output_size)

        self.fc_input_size = self.output_size * out_channels_list[-1]
        self.fc1 = nn.Linear(self.fc_input_size, embedding_size)

        
    def make_layers(self, block, out_channels_list, kernel_sizes, strides):
        layers = []
        num_layers = len(out_channels_list)
        for i in range(0, num_layers):
            layers.append(block(self.in_channels, out_channels_list[i], kernel_sizes[i], stride=strides[i]))
            self.in_channels = out_channels_list[i]
            padding_size = (self.kernel_sizes[i] - 1)/2 
            convolution_stride = 1

            self.output_size = int((self.output_size + 2 * padding_size - kernel_sizes[i]) /convolution_stride + 1) ## convolution layer output_size 
            self.output_size = int((self.output_size  - kernel_sizes[i])/strides[i] + 1) ## maxpool output_size 
        return nn.Sequential(*layers)

    
    def forward(self, x):
        out = x.view(x.size(0), 1,-1) ## reshape numbatch * num_dim to numbatch * num_in_channel * num_dim 
        out = self.layers_block1(out)
        out = out.view(out.size(0), -1)
        # out_size = out.size(0)
        out = self.fc1(out)
        return out
            


class outputLayer(nn.Module):
    def __init__(self, embedding_size, linear_output_size=32, binary_output_size=32):
        super(outputLayer, self).__init__()
        self.linear1 = nn.Linear(embedding_size, linear_output_size) 
        self.linear2 = nn.Linear(embedding_size, binary_output_size)
        self.relu = nn.ReLU()
        self.dense1_bn = nn.BatchNorm1d(1)  
        
    def forward(self, x):
        linear1_out = self.linear1(x)
        print(linear1_out.shape)
        survival_out  = self.dense1_bn(linear1_out[:,0:1])
        linear2_out = self.linear2(x)
        binary_output = self.relu(linear2_out)
        out = torch.cat((survival_out, linear1_out[:,1:], binary_output),1) # define concat
        print(out.shape)
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
        self.dense1_bn = nn.BatchNorm1d(hidden_size)  # batchnorm takes both channels or hidden size as input
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


class Net(nn.Module):
    """
    This is the standard way to define your own network in PyTorch. You typically choose the components
    (e.g. LSTMs, linear layers etc.) of your network in the __init__ function. You then apply these layers
    on the input step-by-step in the forward function. You can use torch.nn.functional to apply functions

    such as F.relu, F.sigmoid, F.softmax, F.max_pool2d. Be careful to ensure your dimensions are correct after each
    step. You are encouraged to have a look at the network in pytorch/nlp/model/net.py to get a better sense of how
    you can go about defining your own network.

    The documentation for all the various components available o you is here: http://pytorch.org/docs/master/nn.html
    """

    def __init__(self, params):
        """
        We define an convolutional network that predicts the sign from an image. The components
        required are:

        - an embedding layer: this layer maps each index in range(params.vocab_size) to a params.embedding_dim vector
        - lstm: applying the LSTM on the sequential input returns an output for each token in the sentence
        - fc: a fully connected layer that converts the LSTM output for each token to a distribution over NER tags

        Args:
            params: (Params) contains num_channels
        """
        super(Net, self).__init__()
        self.num_channels = params.num_channels

        # each of the convolution layers below have the arguments (input_channels, output_channels, filter_size,
        # stride, padding). We also include batch normalisation layers that help stabilise training.
        # For more details on how to use these layers, check out the documentation.
        self.conv1 = nn.Conv2d(3, self.num_channels, 3, stride=1, padding=1)
        self.bn1 = nn.BatchNorm2d(self.num_channels)
        self.conv2 = nn.Conv2d(self.num_channels, self.num_channels * 2, 3, stride=1, padding=1)
        self.bn2 = nn.BatchNorm2d(self.num_channels * 2)
        self.conv3 = nn.Conv2d(self.num_channels * 2, self.num_channels * 4, 3, stride=1, padding=1)
        self.bn3 = nn.BatchNorm2d(self.num_channels * 4)

        # 2 fully connected layers to transform the output of the convolution layers to the final output
        self.fc1 = nn.Linear(8 * 8 * self.num_channels * 4, self.num_channels * 4)
        self.fcbn1 = nn.BatchNorm1d(self.num_channels * 4)
        self.fc2 = nn.Linear(self.num_channels * 4, 6)
        self.dropout_rate = params.dropout_rate

    def forward(self, s):
        """
        This function defines how we use the components of our network to operate on an input batch.

        Args:
            s: (Variable) contains a batch of images, of dimension batch_size x 3 x 64 x 64 .

        Returns:
            out: (Variable) dimension batch_size x 6 with the log probabilities for the labels of each image.

        Note: the dimensions after each step are provided
        """
        #                                                  -> batch_size x 3 x 64 x 64
        # we apply the convolution layers, followed by batch normalisation, maxpool and relu x 3
        s = self.bn1(self.conv1(s))                         # batch_size x num_channels x 64 x 64
        s = F.relu(F.max_pool2d(s, 2))                      # batch_size x num_channels x 32 x 32
        s = self.bn2(self.conv2(s))                         # batch_size x num_channels*2 x 32 x 32
        s = F.relu(F.max_pool2d(s, 2))                      # batch_size x num_channels*2 x 16 x 16
        s = self.bn3(self.conv3(s))                         # batch_size x num_channels*4 x 16 x 16
        s = F.relu(F.max_pool2d(s, 2))                      # batch_size x num_channels*4 x 8 x 8

        # flatten the output for each image
        s = s.view(-1, 8 * 8 * self.num_channels * 4)             # batch_size x 8*8*num_channels*4

        # apply 2 fully connected layers with dropout
        s = F.dropout(F.relu(self.fcbn1(self.fc1(s))),
                      p=self.dropout_rate, training=self.training)    # batch_size x self.num_channels*4
        s = self.fc2(s)                                     # batch_size x 6

        # apply log softmax on each image's output (this is recommended over applying softmax
        # since it is numerically more stable)
        return F.log_softmax(s, dim=1)


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
    torch.isnan(x).any() or x.eq(float('inf')).any() or x.eq(float('-inf')).any()


def isfinite(x):
    """
    Quick pytorch test that there are no nan's or infs.

    note: torch now has torch.isnan
    url: https://gist.github.com/wassname/df8bc03e60f81ff081e1895aabe1f519
    """
    not_inf = ((x + 1) != x)
    not_nan = (x == x)
    return not_inf & not_nan


def negative_log_partial_likelihood(censor, risk, debug=False):
    """Return the negative log-partial likelihood of the prediction
    y_true contains the survival time
    risk is the risk output from the neural network
    censor is the vector of inputs that are censored
    regularization is the regularization constant (not used currently in model)

    Uses the torch backend to perform calculations

    Sorts the surv_time by sorted reverse time
    """

    # calculate negative log likelihood from estimated risk
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


def negative_log_partial_likelihood_loss( risk, censor, debug=False):
    '''
    function to change the order of censor and risk 
    '''

    return negative_log_partial_likelihood(censor, risk, debug)

def c_index(predicted_risk, survival):
    # calculate the concordance index
    survival_time, censor = survival[:, 0], survival[:, 1]
    epsilon = 0.001
    partial_hazard = np.exp(-(predicted_risk + epsilon))
    censor = censor.astype(int)
    ci = concordance_index(survival_time, partial_hazard, censor)
    return ci


# maintain all metrics required in this dictionary- these are used in the training and evaluation loops
metrics = {
    'c_index': c_index,
    # could add more metrics such as accuracy for each token type
}






def calculate_loss(labels, net_outputs, loss_fns):
    '''
    define loss function :
    list of loss function suported are : 
    1. remove NA from the labels 
    2. If there are no observations loss = 0, s.t no derivate.
    Need to test if NA then loss are properly scaled.
    '''

    total_loss = torch.zeros(1)

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

    for i in range(len_fns):
        label, net_output, loss_fn = labels[:,i], net_outputs[:, i], loss_fns[i]
        na_inx = ~np.isnan(labels)
        label, net_output = label[na_inx], net_output[na_inx]
        if(len(label) > 1 ):
            loss_curr = loss_fn(net_output, label) # in case of survival label is censor; censored data assumed to sorted based on patient event time.
        else:
            loss_curr = 0
        total_loss = total_loss + loss_curr 


    return total_loss
