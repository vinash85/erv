'''Transfer learning toy example.
1 - Train a simple convnet on the MNIST dataset the first 5 digits [0..4].
2 - Freeze convolutional layers and fine-tune dense layers
   for the classification of digits [5..9].
Get to 99.8% test accuracy after 5 epochs
for the first five digits classifier
and 99.2% for the last five digits after transfer + fine-tuning.
'''

from __future__ import print_function

import datetime
import keras
from keras.datasets import mnist
from keras.models import Sequential, Model, Input
from keras.layers import Dense, Dropout, Activation, Flatten, BatchNormalization
from keras.layers import Conv1D, MaxPooling1D, Conv2D, MaxPooling2D
from keras import backend as K
import tensorflow as tf
import time
import numpy as np
import argparse
import glob
import json
import os

from keras.callbacks import Callback, ModelCheckpoint, TensorBoard, EarlyStopping, TerminateOnNaN
from keras.models import load_model
from tensorflow.python.lib.io import file_io

import data_generator as gn
from lifelines.utils import concordance_index

now = datetime.datetime.now

batch_size = 128
num_classes = 5
epochs = 2
checkpoint_epochs = 5


train_steps = 80

num_epochs = 200

early_stop = 40

train_batch_size = 256

eval_frequency = 2

# input image dimensions
# number of convolutional filters to use
filters = 8
# size of pooling area for max pooling
pool_size = 2
# convolution kernel size
kernel_size = 3

eval_ci = True

timestr = time.strftime("%Y%m%d_%H%M")


# data_dir = "../data/pancancer_all_immune/"
data_dir = "../data/simulation/"
job_dir_prefix = "../results/model/"

# if os

job_dir = job_dir_prefix + "batch_by_type_normalized_" + timestr


def add2stringlist(prefix, List):
    return [prefix + elem for elem in List]


input_shape = 1076
input_shape = 100

eval_files = validation_files = train_files = add2stringlist(data_dir, ["survival.data.txt", "survival.train.txt"])
# train_files = add2stringlist(data_dir, ["TrainingData.txt", "surv.train.txt"])
# validation_files = add2stringlist(data_dir, ["TestData.txt", "surv.test.txt"])
# eval_files = add2stringlist(data_dir, ["EvalData.txt", "surv.eval.txt"])

# train_files = ["/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/TrainingData.txt", "/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/surv.train.txt"]
# validation_files = ["/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/TestData.txt", "/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/surv.test.txt"]

# eval_files = ["/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/EvalData.txt", "/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/surv.eval.txt"]

# if K.image_data_format() == 'channels_first':
# input_shape = (1, img_rows, img_cols)
# else:
# input_shape = (img_rows, img_cols, 1)


def negative_log_partial_likelihood(censor, risk):
    """Return the negative log-partial likelihood of the prediction
    y_true contains the survival time
    risk is the risk output from the neural network
    censor is the vector of inputs that are censored
    regularization is the regularization constant (not used currently in model)

    Uses the Keras backend to perform calculations

    Sorts the surv_time by sorted reverse time
    """

    # calculate negative log likelihood from estimated risk
    epsilon = 0.00001
    risk = K.reshape(risk, [-1])  # flatten
    hazard_ratio = K.exp(risk)

    # cumsum on sorted surv time accounts for concordance
    log_risk = K.log(tf.cumsum(hazard_ratio) + epsilon)
    log_risk = K.reshape(log_risk, [-1])
    uncensored_likelihood = risk - log_risk

    # apply censor mask: 1 - dead, 0 - censor
    censored_likelihood = uncensored_likelihood * censor
    num_observed_events = K.sum(censor)
    neg_likelihood = - K.sum(censored_likelihood) / \
        tf.cast(num_observed_events, tf.float32)
    print(type(neg_likelihood))

    if (neg_likelihood != neg_likelihood):
        print(neg_likelihood)
        print(censor[np.isnan(censor)])
        print(risk[np.isnan(risk)])
        raise ValueError("nan found")

    return neg_likelihood


class ContinuousEval(Callback):
    """Continuous eval callback to evaluate the checkpoint once
       every so many epochs.
    """

    def __init__(self,
                 kmodel,
                 job_dir,
                 eval_ci,
                 eval_frequency,
                 eval_ci_generator,
                 training_ci_generator):
        self.kmodel = kmodel
        self.job_dir = job_dir
        self.eval_ci = eval_ci
        self.eval_ci_generator = eval_ci_generator
        self.eval_frequency = eval_frequency
        self.training_ci_generator = training_ci_generator

    def concordance_metric(self, survival_time, predicted_risk, censor):
        # calculate the concordance index
        epsilon = 0.001
        partial_hazard = np.exp(-(predicted_risk + epsilon))
        censor = censor.astype(int)
        ci = concordance_index(survival_time, partial_hazard, censor)
        return ci

    def on_epoch_begin(self, epoch, logs={}):
        if epoch > 0 and epoch % self.eval_frequency == 0:

            if self.eval_ci:
                # evaluate CI index for evaluation set
                hazard_features, surv_labels = next(self.eval_ci_generator)

                hazard_predict = self.kmodel.predict(hazard_features)
                ci = self.concordance_metric(
                    surv_labels[:, 0], hazard_predict, surv_labels[:, 1])

                # evaluate CI index for training set
                training_hazard_features, training_surv_labels = next(
                    self.training_ci_generator)

                training_hazard_predict = self.kmodel.predict(
                    training_hazard_features)
                ci_training = self.concordance_metric(
                    training_surv_labels[:, 0], training_hazard_predict, training_surv_labels[:, 1])

                print('\nEvaluation epoch[{}] , Concordance Index Training: {:.2f}, Concordance Index Evaluation:{:.2f}]'.format(
                    epoch, ci_training, ci))
                # write out concordance index to a file for graphing later
                with open(os.path.join(self.job_dir, "concordance.tsv"), "a") as myfile:
                    myfile.write('{}\t{:.2f}\t{:.2f}\n'.format(
                        epoch, ci_training, ci))


def train_model_generator(kmodel, loss_fn):

    kmodel.compile(loss=loss_fn,
                   optimizer=keras.optimizers.Adam(
                       lr=learning_rate, clipvalue=0.5, clipnorm=1.0),
                   metrics=['accuracy'])

    t = now()
    # model.fit(x_train, y_train,
    #           batch_size=batch_size,
    #           epochs=epochs,
    #           verbose=1,
    #           validation_data=(x_test, y_test))

    # check input lengths
    # if (len(train_files) != len(validation_files) or len(train_files) != len(validation_files) or len(validation_files) != len(eval_files)):
    #     raise ValueError("Input file lengths do not match")

    # Read feature files
    train_features = gn.readFile(train_files[0])
    valid_features = gn.readFile(validation_files[0])
    eval_features = gn.readFile(eval_files[0])

    # Read labels
    train_labels = gn.readFile(train_files[1])
    valid_labels = gn.readFile(validation_files[1])
    eval_labels = gn.readFile(eval_files[1])

    train_steps_gen, train_input_size, train_generator = gn.generator_survival(
        train_features, train_labels, shuffle=True, batch_size=train_batch_size)
    valid_steps_gen, valid_input_size, val_generator = gn.generator_survival(
        valid_features, valid_labels, shuffle=True, batch_size=train_batch_size)

    # for ci
    _, _, eval_ci_generator = gn.generator_simple(
        eval_features, eval_labels, batch_size=train_batch_size)
    _, _, training_ci_generator = gn.generator_simple(
        train_features, train_labels, batch_size=train_batch_size)

    # Callbacks

    tblog = TensorBoard(
        log_dir=os.path.join(job_dir, 'logs'),
        histogram_freq=0,
        write_graph=True,
        write_images=False,
        embeddings_freq=0
    )
    tonan = TerminateOnNaN()

    concordance_index_eval = ContinuousEval(
        kmodel,
        job_dir,
        eval_ci,
        eval_frequency,
        eval_ci_generator,
        training_ci_generator)

    kmodel.fit_generator(
        generator=train_generator,
        steps_per_epoch=train_steps,
        epochs=num_epochs,
        validation_data=val_generator,
        validation_steps=10,
        verbose=1,  # for tensorboard visualization
        callbacks=[tblog, tonan, concordance_index_eval])
    print('Training time: %s' % (now() - t))
    # score = kmodel.evaluate(x_test, y_test, verbose=0)
    # print('Test score:', score[0])
    # print('Test accuracy:', score[1])


# feature_layers = [
#     Conv1D(filters, kernel_size,
#            padding='valid',
#            input_shape=input_shape),
#     Activation('relu'),
#     Conv1D(filters, kernel_size),
#     Activation('relu'),
#     MaxPooling1D(pool_size=pool_size),
#     Dropout(0.25),
#     Flatten(),
# ]

# classification_layers = [
#     Dense(128),
#     Activation('relu'),
#     Dropout(0.5),
#     Dense(1),
#     Activation('linear')
# ]

keras.backend.clear_session()

embedder = [
    Dense(64, activation='relu', input_shape=(input_shape,)),
    Dropout(0.1),
    Dense(64, activation='relu'),
    Dropout(0.1)
]

survival_layers = [
    Dense(1, input_shape=(input_shape,))
]
survival_model = Sequential(embedder + survival_layers)

survival_model.summary()

# kmodel = Sequential(classification_layers)

# classification_layers = [
#     # Dense(10, input_shape=(input_shape,)),
#     Dense(1, activation='linear')
# ]
# kmodel = Sequential(feature_layers + classification_layers)

# most stable version
#
keras.backend.clear_session()
inputs = Input(shape=(input_shape,), name='encoder_input')
x = inputs
x = Dense(64, activation='relu')(x)
x = Dropout(0.1)(x)
x = BatchNormalization()(x)
x = Dense(64, activation='relu')(x)
x = Dropout(0.1)(x)
embedding = BatchNormalization()(x)
embedding = Dense(1, activation='linear')(embedding)
embedder = Model(inputs, embedding, name='encoder_output')
embedder.summary()


keras.backend.clear_session()
inputs = Input(shape=(input_shape,), name='encoder_input')
x = inputs
x = Dense(64, activation='relu')(x)
x = Dropout(0.1)(x)
x = BatchNormalization()(x)
x = Dense(64, activation='relu')(x)
x = Dropout(0.1)(x)
embedding = BatchNormalization()(x)
embedder = Model(inputs, embedding, name='encoder_output')
embedder.summary()

# embedder(inputs)

# survival layer
survival_inputs = Input(shape=(64,), name='survival_embedding')
survival_output = Dense(1, activation='linear')(survival_inputs)
survival_layer = Model(survival_inputs, survival_output, name='survival_layer')
survival_layer.summary()

# stack to create survival layer
survival_model = Model(inputs, survival_layer(embedder(inputs)), name="survial_model")
survival_model.summary()
# loss_fn = negative_log_partial_likelihood

# train model for 5-digit classification [0..4]
learning_rate = .0001
num_epochs = 30
train_model_generator(survival_model, loss_fn=negative_log_partial_likelihood)


if False:
    # freeze feature layers and rebuild model
    for l in feature_layers:
        l.trainable = False

    # transfer: train dense layers for new classification task [5..9]
    train_model(model,
                (x_train_gte5, y_train_gte5),
                (x_test_gte5, y_test_gte5), num_classes)

    def layer_output(model, train):
        x_train = train.reshape((train.shape[0],) + input_shape)
        x_train = x_train.astype('float32')
        x_train /= 255
        print('x_train shape:', x_train.shape)
        print(x_train.shape[0], 'train samples')
        return(intermediate_layer_model.predict(x_train))

    from keras.models import Model
    layer_name = 'dense_1'
    intermediate_layer_model = Model(inputs=model.input,
                                     outputs=model.get_layer(layer_name).output)
    last_layer_output = layer_output(intermediate_layer_model, x_test_gte5)

        kmodel.compile(loss=loss_fn,
                       optimizer=keras.optimizers.Adam(
                           lr=learning_rate, clipvalue=0.5, clipnorm=1.0),
                       metrics=['accuracy'])

    last_layer_output = kmodel.predict(aa[0])
    xx = negative_log_partial_likelihood(aa[1], last_layer_output)
    censor, risk = aa[1], last_layer_output

    import keras.backend as K
    import numpy as np

    def get_grads(kmodel, x, y):
        weights = kmodel.trainable_weights  # weight tensors
        gradients = kmodel.optimizer.get_gradients(kmodel.total_loss, weights)  # gradient tensors
        # input_tensors = kkmodel.inputs
        input_tensors = kmodel.inputs + kmodel.sample_weights + kmodel.targets + [K.learning_phase()]
        get_gradients = K.function(inputs=input_tensors, outputs=gradients)
        inputs = [x, np.ones(len(x)), y, 0]
        grads = get_gradients(inputs)
        return(grads)

    aa = next(train_generator)
    out = get_grads(survival_model, aa[0], aa[1])
    # layer_dict[layer_name].get_weights()

    # def get_gradients(model, train):
    #     weights = model.trainable_weights  # weight tensors
    #     gradients = model.optimizer.get_gradients(model.total_loss, weights)
    #     grad_out_model = Model(inputs=model.input,
    #                            outputs=gradients)
    #     return(grad_out_model.predict(train))

    # aa = get_gradients(kmodel, aa[0])
