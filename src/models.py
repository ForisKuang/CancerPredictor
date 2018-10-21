import torch.nn as nn
from torch.nn import functional as F
import torch

class FeedForwardNN(nn.Module):
    def __init__(self, input_len):
        super(FeedForwardNN, self).__init__()
        # Model definition
        self.fc1 = nn.Linear(input_len, 5200)
        self.fc7 = nn.Linear(5200, 2400)
        self.fc2 = nn.Linear(2400, 1200)
        self.fc3 = nn.Linear(1200, 840)
        self.fc4 = nn.Linear(840, 100)
        self.fc5 = nn.Linear(100, 20)
        self.fc6 = nn.Linear(20, 100)

    def forward(self, x):
        x = x.view(-1, self.num_flat_features(x))
        x = F.log_softmax(self.fc1(x), dim=0)
        x = F.log_softmax(self.fc7(x), dim=0)
        x = F.log_softmax(self.fc2(x), dim=0)
        x = F.log_softmax(self.fc3(x), dim=0)
        x = F.log_softmax(self.fc4(x), dim=0)
        x = F.log_softmax(self.fc5(x), dim=0)
        x = self.fc6(x)
        return x

    def num_flat_features(self, x):
        size = x.size()[1:]
        num_features = 1
        for s in size:
            num_features *= s
        return num_features

class LogisticRegression(nn.Module):
    def __init__(self, input_len):
        super(LogisticRegression, self).__init__()
        # Model definition
        self.fc1 = nn.Linear(input_len, 100)

    def forward(self, x):
        x = x.view(-1, self.num_flat_features(x))
        x = self.fc1(x)
        return x

    def num_flat_features(self, x):
        size = x.size()[1:]
        num_features = 1
        for s in size:
            num_features *= s
        return num_features


class RecurrentNN(nn.Module):
    def __init__(self, input_len, hidden_size):
        super(RecurrentNN, self).__init__()
        self.i2h = nn.Linear(input_len + hidden_size, hidden_size)
        self.i2o = nn.Linear(input_len + hidden_size, 100)
        self.o2o = nn.Linear(hidden_size + 100, 100)
        self.dropout = nn.Dropout(0.1)
        self.softmax = nn.LogSoftmax(dim=1)
    
    def forward(self, input, hidden):
        input_combined = torch.cat((input, hidden), 1)
        hidden = self.i2h(input_combined)
        output = self.i2o(input_combined)
        output_combined = torch.cat((hidden, output), 1)
        output = self.o2o(output_combined)
        ouput = self.dropout(output)
        output = self.softmax(output)
        return output, hidden;

    def initHidden(self):
        return torch.zeros(1, self.hidden_size) 

    def num_flat_features(self, x):
        size = x.size()[1:]
        num_features = 1
        for s in size:
            num_features *= s
        return num_features
