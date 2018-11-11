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
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv1D, MaxPooling1D, Conv2D, MaxPooling2D
from keras import backend as K
import tensorflow as tf
import time

import argparse
import glob
import json
import os

from keras.callbacks import Callback, ModelCheckpoint, TensorBoard, EarlyStopping
from keras.models import load_model
from tensorflow.python.lib.io import file_io

import data_generator as gn

now = datetime.datetime.now

batch_size = 128
num_classes = 5
epochs = 2
checkpoint_epochs = 5


train_steps = 80

num_epochs = 10

early_stop = 40

train_batch_size = 256

eval_frequency = 10

input_shape = 1076
# input image dimensions
img_rows, img_cols = 28, 28
# number of convolutional filters to use
filters = 8
# size of pooling area for max pooling
pool_size = 2
# convolution kernel size
kernel_size = 3


timestr = time.strftime("%Y%m%d_%H%M")

job_dir = "/home/as892/project/icb/results/models/batch_by_type_normalized_" + timestr

train_files = ["/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/TrainingData.txt", "/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/surv.train.txt"]
validation_files = ["/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/TestData.txt", "/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/surv.test.txt"]

eval_files = ["/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/EvalData.txt", "/home/pzs2/capstone/proj/TCGA_processed/pancancer_all_immune/surv.eval.txt"]
# if K.image_data_format() == 'channels_first':
# input_shape = (1, img_rows, img_cols)
# else:
# input_shape = (img_rows, img_cols, 1)


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
    tblog = TensorBoard(
        log_dir=os.path.join(job_dir, 'logs'),
        histogram_freq=0,
        write_graph=True,
        write_images=False,
        embeddings_freq=0
    )

    kmodel.fit_generator(
        generator=train_generator,
        steps_per_epoch=train_steps,
        epochs=num_epochs,
        validation_data=val_generator,
        validation_steps=10,
        verbose=1,  # for tensorboard visualization
        callbacks=[tblog])
    print('Training time: %s' % (now() - t))
    # score = kmodel.evaluate(x_test, y_test, verbose=0)
    # print('Test score:', score[0])
    # print('Test accuracy:', score[1])


# for gen in train_generator:
    # aa = gen


# define two groups of layers: feature (convolutions) and classification (dense)


# feature_layers = [
#     Conv2D(filters, kernel_size,
#            padding='valid',
#            input_shape=input_shape),
#     Activation('relu'),
#     Conv2D(filters, kernel_size),
#     Activation('relu'),
#     MaxPooling2D(pool_size=pool_size),
#     Dropout(0.25),
#     Flatten(),
# ]

# classification_layers = [
#     Dense(128),
#     Activation('relu'),
#     Dropout(0.5),
#     Dense(num_classes),
#     Activation('softmax')
# ]

# # create complete model
# model = Sequential(feature_layers + classification_layers)


# # train model for 5-digit classification [0..4]
# train_model(model,
#             (x_train_lt5, y_train_lt5),
    # (x_test_lt5, y_test_lt5), num_classes)


# new model


feature_layers = [
    Conv1D(filters, kernel_size,
           padding='valid',
           input_shape=input_shape),
    Activation('relu'),
    Conv1D(filters, kernel_size),
    Activation('relu'),
    MaxPooling1D(pool_size=pool_size),
    Dropout(0.25),
    Flatten(),
]

classification_layers = [
    Dense(128),
    Activation('relu'),
    Dropout(0.5),
    Dense(1),
    Activation('linear')
]


feature_layers = [
    Dense(128, activation='relu', input_shape=(input_shape,)),
    # Dropout(0.5),
    Dense(128, activation='relu')
    # Dropout(0.5)


]

classification_layers = [
    Dense(1, input_shape=(input_shape,)),
    Activation('relu')
]


kmodel = Sequential(classification_layers)

classification_layers = [
    # Dense(10, input_shape=(input_shape,)),
    Dense(1, activation='linear')
]
kmodel = Sequential(feature_layers + classification_layers)

# most stable version
keras.backend.clear_session()
inputs = Input(shape=(input_shape,), name='encoder_input')
x = inputs
x = Dense(64, activation='relu')(x)
x = Dropout(0.1)(x)
x = Dense(64, activation='relu')(x)
x = Dropout(0.1)(x)
preds = Dense(1, activation='linear')(x)

kmodel = Model(inputs, preds, name='encoder_output')
kmodel.summary()
loss_fn = gn.negative_log_partial_likelihood
# loss_fn = keras.losses.mean_squared_error

# train model for 5-digit classification [0..4]
learning_rate = .0001

train_model_generator(kmodel, loss_fn)


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


out = get_grads(kmodel, aa[0], aa[1])

aa = get_gradients(kmodel)


def get_gradients(model, train):
    weights = model.trainable_weights  # weight tensors
    gradients = model.optimizer.get_gradients(model.total_loss, weights)
    grad_out_model = Model(inputs=model.input,
                           outputs=gradients)
    return(grad_out_model.predict(train))


aa = get_gradients(kmodel, aa[0])
