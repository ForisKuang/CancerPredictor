from gdc_lib import download_files, fetch_patient_data
import utils
import ml_utils
import os
import sys
import numpy as np
import torch.nn as nn
import torch
import torch.utils.data as data_utils
from torch import optim
from models import *

download_files_filter = {
                            "op": "and",
                            "content": [
                                {
                                    "op": "in",
                                    "content": {
                                        "field": "files.access",
                                        "value": ["open"]
                                    }
                                },
                                {
                                    "op": "in",
                                    "content": {
                                        "field": "files.data_category",
                                        "value": ["Simple Nucleotide Variation"]
                                    }
                                },
                                {
                                    "op": "in",
                                    "content": {
                                        "field": "files.data_type",
                                        "value": ["Masked Somatic Mutation"]
                                    }
                                },
                                {
                                    "op": "in",
                                    "content": {
                                        "field": "files.cases.project.program.name",
                                        "value": ["TCGA"]
                                    }
                                }
                            ]
                        }

patient_data_fields = [
                            'case_id',
                            'project.program.name',
                            'project.project_id',
                            'primary_site',
                            'diagnoses.age_at_diagnosis',
                            'diagnoses.tumor_stage',
                            'demographic.gender',
                            'demographic.race',
                            'exposures.alcohol_history',
                            'exposures.cigarettes_per_day',
                            'exposures.years_smoked'
                        ]

patient_data_filter = {
                            "op": "in",
                            "content": {
                                "field": "project.program.name",
                                "value": ["TCGA"]
                            }
                      }

cuda = torch.device('cuda')

def main():
    download_files_flag = raw_input("Do you want to download files? (Y/N) ")
    if (download_files_flag.upper().startswith("Y")):
        num_files = int(raw_input("Number of files you would like to download "))
        download_files(download_files_filter, num_files)
    tar_files_flag = raw_input("Do you want to tar files? (Y/N) ")
    if (tar_files_flag.upper().startswith("Y")):
        utils.tar_files("datasets")

    write_data_to_file_flag = raw_input("Do you want to write the data to a file? (Y/N) ")
    if (write_data_to_file_flag.upper().startswith("Y")):
        patient_data = fetch_patient_data(patient_data_fields, 5000, data_format="TSV", filters=patient_data_filter)
        print ("Got patient data")
        case_ids = utils.get_case_ids("datasets", patient_data)
        print("Case_ids finished")
        eighty_percent_of_case_ids = round(0.8*len(case_ids), 0)
        train_case_ids = {}
        test_case_ids = {}
        count = 0
        for k, v in case_ids.iteritems():
            if (count < eighty_percent_of_case_ids):
                train_case_ids[k] = v
            else:
                test_case_ids[k] = v
            count += 1

        print("Creating train and test case_ids finished")
        indices = utils.get_useful_indices_list("gen_seq/gene_annotation", "gen_seq/Oncogenes.csv")
        print("Indices Finished")

        Xdata, Ydata, Xtest, Ytest = ml_utils.construct_data(indices, train_case_ids, test_case_ids)

        utils.write_data_to_files(Xdata, Ydata, Xtest, Ytest)
    
    
    Xdata, Ydata, Xtest, Ytest = utils.get_numpy_data_from_files("datasets/Xdata.csv", \
        "datasets/Ydata.csv", "datasets/Xtest.csv", "datasets/Ytest.csv")

    print("Xdata shape:", Xdata.shape)
    print("Ydata shape:", Ydata.shape)
    print("Xtest shape:", Xtest.shape)
    print("Ytest shape:", Ytest.shape)

    tensor_train_x = torch.from_numpy(Xdata)
    tensor_train_y = torch.from_numpy(Ydata)

    tensor_test_x = torch.from_numpy(Xtest)
    tensor_test_y = torch.from_numpy(Ytest)

    train_dataset = data_utils.TensorDataset(tensor_train_x, tensor_train_y)
    test_dataset = data_utils.TensorDataset(tensor_test_x, tensor_test_y)


    models_decision = raw_input("Which model do you want to try? (logistic, ff, rnn) ")
    
    net = None
    criterion = None
    optimizer = None
    epochs = 1000
    train_dataloader = data_utils.DataLoader(train_dataset, batch_size=100, shuffle=True, num_workers=2)
    test_dataloader = data_utils.DataLoader(test_dataset, batch_size=1, shuffle=False, num_workers=2)

    if (models_decision.lower().startswith('l')):
        net = LogisticRegression(17335)
        net = net.to(cuda)

        criterion = nn.CrossEntropyLoss()
        optimizer = optim.SGD(net.parameters(), lr=0.01, momentum=0.9, weight_decay=5e-4)
        for epoch in range(epochs):
            print("On epoch " + str(epoch) + " out of " + str(epochs))
            ml_utils.train(net, train_dataloader, optimizer, criterion, epoch)
            print("Testing train")
            ml_utils.test(net, train_dataloader)
            print("Testing test")
            ml_utils.test(net, test_dataloader)
    
    elif (models_decision.lower().startswith('f')):
        net = FeedForwardNN(17335)
        net = net.to(cuda)

        criterion = nn.CrossEntropyLoss()
        optimizer = optim.SGD(net.parameters(), lr=0.01, momentum=0.9, weight_decay=5e-4)
        for epoch in range(epochs):
            print("On epoch " + str(epoch) + " out of " + str(epochs))
            ml_utils.train(net, train_dataloader, optimizer, criterion, epoch)
            print("Testing train")
            ml_utils.test(net, train_dataloader)
            print("Testing test")
            ml_utils.test(net, test_dataloader)
    
    elif (models_decision.lower().startswith('r')):
        net = RecurrentNN(17335, 5200)
        net = net.to(cuda)

        criterion = nn.CrossEntropyLoss()
        optimizer = optim.SGD(net.parameters(), lr=0.01, momentum=0.9, weight_decay=5e-4)
    else:
        print("Not a valid input")
        raise SystemExit


    #w = ml_utils.closed_form_linreg(Xdata, Ydata, 0.001)
    #print("Train Class Error is " + str(ml_utils.compute_class_error(Ydata, np.dot(Xdata, w))))
    #print("Test Class Error is " + str(ml_utils.compute_class_error(Ytest, np.dot(Xtest, w))))
    #print("Sq Error is " + str(ml_utils.compute_Sq_error(Ydata, Xdata, w)))


if __name__ == '__main__':
    main()
