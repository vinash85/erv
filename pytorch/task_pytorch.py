import torch
import torch.nn as nn
import torchvision
import torchvision.transforms as transforms
# from ... import data_generator as gn
import data_generator as gn


# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
cudnn.benchmark = True

# Hyper parameters
num_epochs = 5
num_classes = 1
batch_size = 100
learning_rate = 0.001


# read features


# Read feature files

# data_dir = "/home/as892/project/icb/data/"
data_dir = "../../data/ssgsea/"
job_dir_prefix = "../../results/model/"


job_dir = job_dir_prefix + "batch_ssgsea_" + timestr


def add2stringlist(prefix, List):
    return [prefix + elem for elem in List]


train_files = add2stringlist(data_dir, ["tcga_ssgsea_train.txt", "tcga_survival_train.txt"])
validation_files = add2stringlist(data_dir, ["tcga_ssgsea_test.txt", "tcga_survival_test.txt"])
eval_files = add2stringlist(data_dir, ["tcga_ssgsea_eval.txt", "tcga_survival_eval.txt"])


train_features = gn.readFile(train_files[0])
test_features = gn.readFile(test_files[0])

# Read labels
train_labels = gn.readFile(train_files[1])
test_labels = gn.readFile(test_files[1])

train_steps_gen, train_input_size, train_generator = gn.generator_survival(
    train_features, train_labels, shuffle=True, batch_size=train_batch_size)
test_steps_gen, test_input_size, val_generator = gn.generator_survival(
    test_features, test_labels, shuffle=True, batch_size=train_batch_size)


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


model = NeuralNet(input_size, hidden_size, num_classes).to(device)

# model = ConvNet(num_classes).to(device)

# Loss and optimizer
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Train the model
total_step = len(train_loader)
for epoch in range(num_epochs):
    for i, (images, labels) in enumerate(train_loader):
        images = images.to(device)
        labels = labels.to(device)

        # Forward pass
        outputs = model(images)
        loss = criterion(outputs, labels)

        # Backward and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (i + 1) % 100 == 0:
            print ('Epoch [{}/{}], Step [{}/{}], Loss: {:.4f}'
                   .format(epoch + 1, num_epochs, i + 1, total_step, loss.item()))

# Test the model
model.eval()  # eval mode (batchnorm uses moving mean/variance instead of mini-batch mean/variance)
with torch.no_grad():
    correct = 0
    total = 0
    for images, labels in test_loader:
        images = images.to(device)
        labels = labels.to(device)
        outputs = model(images)
        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

    print('Test Accuracy of the model on the 10000 test images: {} %'.format(100 * correct / total))

# Save the model checkpoint
torch.save(model.state_dict(), 'model.ckpt')
