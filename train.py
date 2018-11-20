import torch
import torch.nn as nn
import torchvision
import torchvision.transforms as transforms
# from ... import data_generator as gn
import data_generator_pytorch as gn
import datetime
import time
import metrics

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
# cudnn.benchmark = True

# Hyper parameters
num_epochs = 200
num_classes = 1
batch_size = 100
learning_rate = 0.001
hidden_size = 256
now = datetime.datetime.now
t = now()
timestr = time.strftime("%Y%m%d_%H%M")

train_batch_size = 256
# read features


# Read feature files


# define concordance index


# define the loss function


def negative_log_partial_likelihood(censor, risk):
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
    risk = torch.reshape(risk, [-1])  # flatten
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
    # print(type(neg_likelihood))

    if (neg_likelihood != neg_likelihood):
        print(neg_likelihood)
        print(censor[np.isnan(censor)])
        print(risk[np.isnan(risk)])
        raise ValueError("nan found")

    return neg_likelihood


# Fully connected neural network with one hidden layer


class NeuralNet(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(NeuralNet, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
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
    with tqdm(total=len(dataloader)) as t:
        for i, (features, survival) in zip(range(total_step), train_generator):
            labels_batch = survival[:, 1]
            # move to GPU if available
            if params.cuda:
                train_batch, labels_batch = train_batch.cuda(async=True), labels_batch.cuda(async=True)
            # convert to torch Variables
            train_batch, labels_batch = Variable(train_batch), Variable(labels_batch)

            # compute model output and loss
            output_batch = model(train_batch)
            loss = loss_fn(output_batch, labels_batch)

            # clear previous gradients, compute gradients of all variables wrt loss
            optimizer.zero_grad()
            loss.backward()

            # performs updates using calculated gradients
            optimizer.step()

            # Evaluate summaries only once in a while
            if i % params.save_summary_steps == 0:
                # extract data from torch Variable, move to cpu, convert to numpy arrays
                output_batch = output_batch.data.cpu().numpy()
                labels_batch = labels_batch.data.cpu().numpy()
                c_index = metrics.concordance_metric(output_batch, survival)

                # compute all metrics on this batch
                summary_batch = {'c_index': c_index}
                summary_batch['loss'] = loss.data[0]
                summ.append(summary_batch)

            # update the average loss
            loss_avg.update(loss.data[0])

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

    best_val_acc = 0.5  # for confidence

    for epoch in range(params.num_epochs):
        # Run one epoch
        logging.info("Epoch {}/{}".format(epoch + 1, params.num_epochs))

        # compute number of batches in one epoch (one full pass over the training set)
        train(model, optimizer, loss_fn, train_dataloader, metrics, params)

        # Evaluate for one epoch on validation set
        val_metrics = evaluate(model, loss_fn, val_dataloader, metrics, params)

        val_acc = val_metrics['c_index']
        is_best = val_acc >= best_val_acc

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

    # use GPU if available
    params.cuda = torch.cuda.is_available()

    # Set the random seed for reproducible experiments
    torch.manual_seed(230)
    if params.cuda:
        torch.cuda.manual_seed(230)

    # Set the logger
    utils.set_logger(os.path.join(args.model_dir, 'train.log'))

    # Create the input data pipeline
    logging.info("Loading the datasets...")

    # # data_dir = "/home/as892/project/icb/data/"
    # data_dir = "../../data/ssgsea/"
    # # data_dir = "../../results/simulation/"
    # job_dir_prefix = "../../results/model/"

    # job_dir = job_dir_prefix + "pytorch_ssgsea_" + timestr

    def add2stringlist(prefix, List):
        return [prefix + elem for elem in List]

    # train_files = add2stringlist(data_dir, ["survival.data.txt", "survival.train.txt"])
    # test_files = add2stringlist(data_dir, ["survival.data.txt", "survival.train.txt"])
    train_files = add2stringlist(data_dir, ["tcga_ssgsea_train.txt", "tcga_survival_train.txt"])
    test_files = add2stringlist(data_dir, ["tcga_ssgsea_test.txt", "tcga_survival_test.txt"])

    train_features = gn.readFile(train_files[0])
    test_features = gn.readFile(test_files[0])

    # Read labels
    train_labels = gn.readFile(train_files[1])
    test_labels = gn.readFile(test_files[1])

    train_steps_gen, train_input_size, train_generator = gn.generator_survival(
        train_features, train_labels, shuffle=True, batch_size=train_batch_size)
    test_steps_gen, test_input_size, val_generator = gn.generator_survival(
        test_features, test_labels, shuffle=True, batch_size=train_batch_size)

    input_size = train_input_size

    # fetch dataloaders
    # dataloaders = data_loader.fetch_dataloader(['train', 'val'], args.data_dir, params)
    train_dl = train_generator
    val_dl = val_generator

    logging.info("- done.")

    # Define the model and optimizer
    model = net.FCN(params).cuda() if params.cuda else net.FCN(params)
    optimizer = optim.Adam(model.parameters(), lr=params.learning_rate)

    # fetch loss function and metrics
    loss_fn = net.negative_log_partial_likelihood
    metrics = net.metrics

    # Train the model
    logging.info("Starting training for {} epoch(s)".format(params.num_epochs))
    train_and_evaluate(model, train_dl, val_dl, optimizer, loss_fn, metrics, params, args.model_dir,
                       args.restore_file)


# model = NeuralNet(input_size, hidden_size, num_classes).double().to(device)

# model = ConvNet(num_classes).to(device)

# Loss and optimizer
# criterion = nn.CrossEntropyLoss()
# optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Train the model
# total_step = len(train_loader)
total_step = train_steps_gen
for epoch in range(num_epochs):
    for i, (features, survival) in zip(range(total_step), train_generator):
        features = (torch.from_numpy(features)).to(device)
        censor = (torch.from_numpy(survival[:, 1])).double().to(device)

        # Forward pass
        outputs = model(features)
        # loss = criterion(outputs, labels)
        loss = negative_log_partial_likelihood(censor, outputs)

        # Backward and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (i + 1) % 1 == 0:
            c_index = metrics.concordance_metric(outputs.detach().cpu().numpy(), survival)
            print ('Epoch [{}/{}], Step [{}/{}], Loss: {:.4f}, c_index: {:.4f} '
                   .format(epoch + 1, num_epochs, i + 1, total_step, loss.item(), c_index))

# Test the model
model.eval()  # eval mode (batchnorm uses moving mean/variance instead of mini-batch mean/variance)
with torch.no_grad():
    correct = 0
    total = 0
    features = (torch.from_numpy(test_features)).to(device)
    # censor = (torch.from_numpy(test_labels[:, 1])).to(device)
    outputs = model(features)
    c_index = metrics.concordance_metric(outputs.cpu().detach().numpy(), test_labels)
    # for features, survival in train_generator:
    #     features = (torch.from_numpy(features)).to(device)
    #     censor = (torch.from_numpy(survival[:, 1])).to(device)
    #     outputs = model(features)
    #     c_index = metrics.concordance_metric(outputs.detach().numpy(), survival)

    print('Test c-index of the model on the test data: {} %'.format(c_index))

# Save the model checkpoint
torch.save(model.state_dict(), 'model.ckpt')
