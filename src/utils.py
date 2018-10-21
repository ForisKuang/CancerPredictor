import os, gzip, csv
import tarfile
from subprocess import call
from gdc_lib import fetch_patient_data
from numpy import genfromtxt
import numpy as np

def get_file_ids(directory):
    """
        params: directory
                    the directory in which to traverse to begin getting file ids
        returns: file_ids
                    a map whose key values are unique ids of the files and the
                    value is a tuple with the first value being the project id
                    and the second value being the scan type
    """
    file_ids = {}
    for filename in os.listdir(directory):
        split_filename = filename.split('.')
        file_ids[split_filename[3]] = (split_filename[1], split_filename[2])
    return file_ids

def get_case_ids(directory, patient_data):
    """
        params: directory
                    the directory in which to traverse to aggregate the case ids
        returns: case_ids
                    a map of case ids to ({(gene, indices)*}, age, stage)
    """
    case_ids = {}
    ids_from_patient_data = get_case_ids_from_patient_data(patient_data)
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if ("maf.gz" in filename):
                split_filename = filename.split('.')
                ids_from_file = get_case_ids_from_file(directory + '/' + split_filename[3] + '/' + filename)
                for key in ids_from_patient_data.keys():
                    if (key in ids_from_file.keys()):
                        case_ids[key] = (ids_from_file[key], ids_from_patient_data[key][0], ids_from_patient_data[key][1])
    return case_ids

def get_case_ids_from_patient_data(patient_data):
    """
        params: content
                    Response content from fetching data
        returns: case_ids
                    A map of case ids to a tuple of (age, tumor_stage)
                
    """
    first_line = True
    case_id_idx = None
    age_at_diagnosis_idx = None
    tumor_stage_idx = None
    case_ids = {}
    split_data = patient_data.splitlines()
    for line in split_data:
        split_line = line.split('\t')
        if (first_line):
            for i in range(len(split_line)):
                if ("case_id" in split_line[i]):
                    case_id_idx = i
                elif ("tumor_stage" in split_line[i]):
                    tumor_stage_idx = i
                elif ("age_at_diagnosis" in split_line[i]):
                    age_at_diagnosis_idx = i
            first_line = False
        else:
            if len(split_line[age_at_diagnosis_idx]) and split_line[tumor_stage_idx] != u'not reported': case_ids[split_line[case_id_idx]] = (int(float(split_line[age_at_diagnosis_idx])), split_line[tumor_stage_idx].replace("stage ", ""))
    return case_ids

def get_case_ids_from_file(filename):
    """
        params: filename
                    File should be of type maf.gz in order to have consistent data
        returns: case_ids
                    A map of case ids to a set of (Gene and indices)
                
    """
    # Read file contents
    with gzip.open(filename) as f:
        stuff = f.read()

        # Parse contents to find correct data
        stuff = stuff.decode()
        stuff = stuff[stuff.find("Hugo"):]
        stuff = stuff.splitlines()[1:]

        # Map in the form of case_id->{(hugo_sym, (start_mut, end_mut))}
        # for non-silent mutations
        case_id_map = {}
        hugo_sym_ind = 0
        start_mut_ind = 5
        end_mut_ind = 6
        mut_type_ind = 8
        case_id_ind = 115
        for line in stuff:
            line = line.split('\t')
            if (line[mut_type_ind] != "Silent"):
                case_id = line[case_id_ind]
                if (case_id not in case_id_map):
                    case_id_map[case_id] = set()

                hugo_sym = line[hugo_sym_ind]
                start_mut = line[start_mut_ind]
                end_mut = line[end_mut_ind]
                case_id_map[case_id].add((hugo_sym, (start_mut, end_mut)))

        return case_id_map

def write_case_id_to_tsv(case_ids, filename):
    """
        params: case_ids
                    The map of case ids to gene mutations, stage, and age
                filename
                    The resulting file name
        output: writes caseids to the filename in tsv format
    """
                
    with open(filename, 'wb') as f:
        for k, v in case_ids.iteritems():
            tab_separated = str(k) + '\t' + str(v[1]) + '\t' + str(v[2]) + '\t' + str(v[0]) + '\n'
            f.write(tab_separated)
            
def find_most_mutations_in_cases(case_id_map):
    max_len = 0
    max_set = None
    for k, v in case_id_map.iteritems():
        if (len(v[0]) > max_len):
            max_len = len(v[0])
            max_set = v[0]
    return max_len, max_set

def tar_files(directory):
    """
        params: directory
                    The directory in which you want to tar files
        output: untarred files of anything that has a .tar.gz tag
    """
    for filename in os.listdir(directory):
        if "tar.gz" in filename:
            print(filename)
            tar = tarfile.open(directory + "/" + filename)
            tar.extractall(path=directory)
            tar.close() 
            os.remove(directory + "/" + filename)

def ungzip(filename):
    with gzip.open(filename) as f:
        val = f.read()
    return val.decode()

def get_useful_indices_list(indFile, geneFile):
   ###
   with open(geneFile, 'r') as g:
   ###
       with open(indFile, 'r') as f:
            ###
            genes = g.read().split(',')
            ###

            line_sep = f.read().split('\n')

            inds = []
            filters = ["gene", "transcript"]
            for line in line_sep:
                if len(line) > 0 and line[0] != '#':
                    splt = line.split()

                    ###
                    gene = splt[8][1:-2]
                    for i in range(8, len(splt)):
                        if (splt[i] == "gene_name"):
                            gene = splt[i+1][1:-2]
                            break
                    ###

                    ### (gene in genes) part
                    if gene in genes and splt[2] not in filters:
                        inds.append((int(splt[3]), int(splt[4])))

            return sorted(inds, key=lambda x: x[0])

def binary_search(lst, data):
    low = 0
    high = len(lst) - 1
    idx = None
    while low <= high:
        mid = low +  ((high - low)/2)
        if (lst[mid][0] <= data[0] and data[0] <= lst[mid][1]) or (lst[mid][0] <= data[1] and data[1] <= lst[mid][1]):
            idx = mid
            break
        else:
            if data[1] < lst[mid][0]:
                high = mid - 1
            else:
                low = mid + 1
    return idx

def write_to_csv(data, filename):
    np.savetxt(filename, data, delimiter=',')
    """
    with open(filename, "w+") as f:
        csvWriter = csv.writer(f)
        csvWriter.writerows(data)
    """

def write_data_to_files(Xdata, Ydata, Xtest, Ytest):
    write_to_csv(Xdata, "datasets/Xdata.csv")
    write_to_csv(Ydata, "datasets/Ydata.csv")
    write_to_csv(Xtest, "datasets/Xtest.csv")
    write_to_csv(Ytest, "datasets/Ytest.csv")

def read_from_csv_to_numpy(filename):
    return genfromtxt(filename, delimiter=',')

def get_numpy_data_from_files(XdataFile, YdataFile, XtestFile, YtestFile):
    Xdata = read_from_csv_to_numpy(XdataFile)
    Ydata = read_from_csv_to_numpy(YdataFile)
    Xtest = read_from_csv_to_numpy(XtestFile)
    Ytest = read_from_csv_to_numpy(YtestFile)

    return Xdata, Ydata, Xtest, Ytest
    
