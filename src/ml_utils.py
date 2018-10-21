import numpy as np
import utils
import torch
import torch.nn as nn
from torch.nn import functional as F
from torch.autograd import Variable
from torch import optim
import math, random

stages = {"i" : 0,
          "ia" : 1,
          "ib": 2,
          "ic": 3,
          "is": 4,
          "ii": 5,
          "iia": 6,
          "iib": 7,
          "iic": 8,
          "iii": 9,
          "iiia": 10,
          "iiib": 11,
          "iiic": 12,
          "iv": 13,
          "iva": 14,
          "ivb": 15,
          "ivc": 16,
          "Error": 17}

cuda = torch.device('cuda')


def feature_map(indices, case):
    """
        Non binary feature map for lin reg
    """
    phi_x = np.asarray([0]*(len(indices) + 1))
    for gene, gene_idx in case[0]:
        idx = utils.binary_search(indices, gene_idx)
        if idx is None:
            continue
        if (gene_idx[0] >= indices[idx][0] and gene_idx[1] <= indices[idx][1]):
            phi_x[i] += (gene_idx[1] - gene_idx[0]) + 1
        elif ((gene_idx[0] >= indices[idx][0] and gene_idx[0] <= indices[idx][1]) \
                or (gene_idx[1] >= indices[idx][0] and gene_idx[1] <= indices[idx][1])):
            if (indices[idx][1] - gene_idx[1] >= 0):
                phi_x[idx] += gene_idx[1] - indices[idx][0]
            else:
                phi_x[idx] += indices[idx][1] - gene_idx[0]
    stages_lst = [0]*18
    if case[2] not in stages:
        stages_lst[17] = 1
        phi_x = np.concatenate((phi_x, stages_lst), axis=0)
    else:
        stages_lst[stages[case[2]]] = 1
        phi_x = np.concatenate((phi_x, stages_lst), axis=0)
    return phi_x

def construct_data(indices, train_case_ids, test_case_ids):
    Xtrain = []
    Ytrain = []
    Xtest = []
    Ytest = []
    for k, v in train_case_ids.iteritems():
        Xtrain.append(feature_map(indices, v))
        val = v[1]//365
        if (val >= 100):
            val = 100
        Ytrain.append(val)
    for k, v in test_case_ids.iteritems():
        Xtest.append(feature_map(indices, v))
        val = v[1]//365
        if (val >= 100):
            val = 100
        Ytest.append(val)
    return np.asarray(Xtrain), np.asarray(Ytrain), np.asarray(Xtest), np.asarray(Ytest)

def closed_form_linreg(Xdata, Ydata, lambda_val):
    """
        params: Xdata
                    an N x d where N is number of data points and d is number of dimensions per data point
                Ydata
                    an N x 1 array with the correct output (ages)
                lambda_val
                    Regularlizer to lower later on
        returns: weight_vector
    """
    N, d = Xdata.shape
    multiply_X = np.dot(Xdata.T, Xdata)
    div_by_N = np.divide(multiply_X, float(N))
    reg = np.multiply(lambda_val, np.identity(d))
    inv_X = np.linalg.inv(div_by_N + reg)
    y_mod = np.divide(np.dot(Xdata.T, Ydata), float(N))
    w = np.dot(inv_X, y_mod)
    return w

def compute_class_error(Y, Yhat):
    errors = 0
    for i in range(len(Yhat)):
        if (abs(round(Yhat[i], 0) - Y[i]) >= 5):
            errors += 1
    print(errors)
    return np.sum(100*errors/len(Yhat))

def compute_sq_error(Y, X, w):
    return np.sum(np.square(Y - np.dot(X, w))/(2*len(X)))
    
def train(net, dataloader, optimizer, criterion, epoch):
    
    total_loss = 0.0
    running_loss = 0.0

    for i, data in enumerate(dataloader, 0):
        inputs, labels = data
        inputs = inputs.type(torch.FloatTensor)
        inputs = F.normalize(inputs)
        inputs = inputs.to(cuda)
        labels = labels.type(torch.LongTensor)
        labels = labels.to(cuda)

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()
        running_loss += loss.item()
        if (i + 1) % 2000 == 0:
            print(running_loss/2000)
            running_loss = 0.0
    print("The total loss is " + str(total_loss/i))
        

def test(net, dataloader):
    correct = 0
    total = 0
    age_count = [0] * 4
    age_correct = [0] * 4
    dataTestLoader = dataloader
    with torch.no_grad():
        for data in dataTestLoader:
            inputs, labels = data
            inputs = inputs.type(torch.FloatTensor)
            inputs = F.normalize(inputs)
            inputs = inputs.to('cuda')
            labels = labels.type(torch.LongTensor)
            labels = labels.to('cuda')
            outputs = net(inputs)
            values, predicted = torch.max(outputs.data, 1)
            total += labels.size(0)
            for i in range(labels.size(0)):
                x = predicted[i][0].item()
                y = labels[i][0].item()
                if (y < 25):
                    age_count[0] += 1
                    if(abs(x - y) <= 5):
                        correct += 1
                        age_correct[0] += 1
                elif (y < 50):
                    age_count[1] += 1
                    if(abs(x - y) <= 5):
                        correct += 1
                        age_correct[1] += 1
                elif (y < 75):
                    age_count[2] += 1
                    if(abs(x - y) <= 5):
                        correct += 1
                        age_correct[2] += 1
                else:
                    age_count[3] += 1
                    if(abs(x - y) <= 5):
                        correct += 1
                        age_correct[3] += 1
                    
        print("Test Accuracy is " + str(float(correct)/float(total)))
        print("Accuracy for 0-24: " + str(float(age_correct[0])/float(age_count[0])))
        print("Accuracy for 25-49: " + str(float(age_correct[1])/float(age_count[1])))
        print("Accuracy for 50-74: " + str(float(age_correct[2])/float(age_count[2])))
        print("Accuracy for 75-100+: " + str(float(age_correct[3])/float(age_count[3])))
