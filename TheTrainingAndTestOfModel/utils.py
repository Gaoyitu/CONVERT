import numpy as np
import pandas as pd

def get_data_and_label(data,label,num_of_labels):

    data_list = data.split(",")
    label_array = np.zeros(shape=(num_of_labels),dtype=np.float32)
    data_array = np.zeros(shape=(data_list.__len__(),21),dtype=np.float32)
    label_array[int(label)] = 1.0
    label_array.resize(1,num_of_labels)
    for i in range(data_list.__len__()):
        data_array[i][int(data_list[i])] = 1.0
    data_array.resize([1,data_array.shape[0],data_array.shape[1]])
    return data_array,label_array

def get_data(data):

    data_list = data.split(",")
    data_array = np.zeros(shape=(data_list.__len__(),21),dtype=np.float32)
    for i in range(data_list.__len__()):
        data_array[i][int(data_list[i])] = 1.0
    data_array.resize([1,data_array.shape[0],data_array.shape[1]])
    return data_array

def get_protein_vector(data):
    data_list = data.split(",")
    data_array = np.zeros(shape=(1,data_list.__len__()),dtype=np.float32)
    for i in range(data_list.__len__()):
        data_array[0][i] = np.float32(data_list[i])
    return data_array
def get_batch_data_and_label(data_list,label_list,num_of_labels,batch_size):

    lengthest_size = 0
    for i in range(len(data_list)):
        if len(data_list[i])>lengthest_size:
            lengthest_size = len(data_list[i])
    lengthest_size = int((lengthest_size+1)/2)

    batch_data_array = np.zeros(shape=(batch_size,lengthest_size,21),dtype=np.float32)
    batch_label_array = np.zeros(shape=(batch_size,num_of_labels),dtype=np.float32)
    for i in range(len(data_list)):
        index = data_list[i].split(",")
        for j in range(len(index)):
            batch_data_array[i][j][int(index[j])] = 1.0
    for i in range(len(label_list)):
        batch_label_array[i][int(label_list[i])] = 1.0
    return batch_data_array,batch_label_array

def get_shuffle_data(file_name):
    f = open(file_name)
    data_list = []
    rep_list = []
    while (True):
        id = f.readline()
        if id == '':
            f.close()
            break
        data = f.readline()
        rep = f.readline()
        data_list.append(''.join(data).strip('\n'))
        rep_list.append(''.join(rep).strip('\n'))
    dictionary = {"data": data_list, "label": rep_list}
    dictionary = pd.DataFrame(dictionary)
    dictionary = dictionary.sample(frac=1)
    data_list = list(dictionary["data"])
    label_list = list(dictionary["label"])
    return data_list,label_list
def over_sampling(data):
    from collections import Counter
    import random
    f = open("data/selectionData_add_rep/train/"+data+".txt")
    write = open("data/experimentalData/train/"+data+".txt", "w")
    label_list = []
    data_list = []
    rep_list = []
    while (True):
        label = f.readline()
        if label == '':
            f.close()
            break
        data = f.readline()
        rep = f.readline()
        label_list.append(label)
        data_list.append(data)
        rep_list.append(rep)
    most_seq_num = Counter(label_list).most_common(1)[0][1]
    print(most_seq_num)
    label_num = "a.1.1.1\n"
    family_seq = []
    family_id = []
    family_rep = []
    one_family_seq = []
    one_family_id = []
    one_family_rep = []
    for i in range(len(label_list)):
        if label_num != label_list[i]:
            label_num = label_list[i]
            family_seq.append(one_family_seq)
            family_id.append(one_family_id)
            family_rep.append(one_family_rep)
            one_family_seq = []
            one_family_id = []
            one_family_rep = []
            one_family_seq.append(data_list[i])
            one_family_id.append(label_list[i])
            one_family_rep.append(rep_list[i])
            if i == len(label_list) - 1:
                family_seq.append(one_family_seq)
                family_id.append(one_family_id)
                family_rep.append(one_family_rep)
        else:
            one_family_seq.append(data_list[i])
            one_family_id.append(label_list[i])
            one_family_rep.append(rep_list[i])
            if i == len(label_list) - 1:
                family_seq.append(one_family_seq)
                family_id.append(one_family_id)
                family_rep.append(one_family_rep)

    for i in range(len(family_seq)):
        one_family_seq = family_seq[i]
        one_family_id = family_id[i]
        one_family_rep = family_rep[i]
        while (len(one_family_seq) < most_seq_num):
            one_family_seq.append(family_seq[i][random.randint(0, len(family_seq[i]) - 1)])
            one_family_id.append(family_id[i][0])
            one_family_rep.append(family_rep[i][0])
        for j in range(len(one_family_rep)):
            write.write(one_family_id[j])
            write.write(one_family_seq[j])
            write.write(one_family_rep[j])
    write.close()

def data_selection(data):
    import random
    f = open("data/sourceData/"+data+".txt")
    f1 = open("data/sourceRepProData/"+data+".txt")
    write_train = open("data/selectionData/train/"+data+".txt", "w")
    write_test = open("data/selectionData/test/"+data+".txt", "w")
    write_experimental = open("data/experimentalRepPro/"+data+".txt","w")
    label_list = []
    data_list = []
    while (True):
        label = f.readline()
        if label == '':
            f.close()
            break
        data = f.readline()
        label_list.append(label)
        data_list.append(data)
    rep_label_list = []
    rep_data_list = []
    while (True):
        label = f1.readline()
        if label == '':
            f1.close()
            break
        data = f1.readline()
        rep_label_list.append(label)
        rep_data_list.append(data)
    label_num = "a.1.1.1\n"
    family_seq = []
    family_id = []
    one_family_seq = []
    one_family_id = []
    for i in range(len(label_list)):
        if label_num != label_list[i]:
            label_num = label_list[i]
            family_seq.append(one_family_seq)
            family_id.append(one_family_id)
            one_family_seq = []
            one_family_id = []
            one_family_seq.append(data_list[i])
            one_family_id.append(label_list[i])
            if i == len(label_list) - 1:
                family_seq.append(one_family_seq)
                family_id.append(one_family_id)
        else:
            one_family_seq.append(data_list[i])
            one_family_id.append(label_list[i])
            if i == len(label_list) - 1:
                family_seq.append(one_family_seq)
                family_id.append(one_family_id)
    for i in range(len(family_id)):
        if len(family_id[i]) > 4:
            rand = random.randint(0, len(family_id[i]) - 1)
            write_experimental.write(rep_label_list[i])
            write_experimental.write(rep_data_list[i])
            for j in range(len(family_id[i])):
                if j == rand:
                    write_test.write(family_id[i][j])
                    write_test.write(family_seq[i][j])

                else:
                    write_train.write(family_id[i][j])
                    write_train.write(family_seq[i][j])
                    write_train.write(rep_data_list[i])
    write_train.close()
    write_test.close()
    write_experimental.close()

def get_test_data(file_name):
    f = open(file_name)
    label_list = []
    data_list = []
    while (True):
        label = f.readline()
        if label == '':
            f.close()
            break
        data = f.readline()
        label_list.append(label)
        data_list.append(data)
    return data_list,label_list

def get_trian_data(train_data_file):
    f = open(train_data_file)
    seq_list = []
    tar_list = []
    while True:
        id = f.readline()
        if id == "":
            f.close()
            break
        x = f.readline()
        y = f.readline()
        seq_list.append(x)
        tar_list.append(y)
    return seq_list,tar_list

def get_batch_datas_and_labels(data_list,label_list,batch_size):

    data_lengthest_size = 0
    label_lengthest_size = 0
    for i in range(len(data_list)):
        if len(data_list[i])>data_lengthest_size:
            data_lengthest_size = len(data_list[i])
    data_lengthest_size = int((data_lengthest_size+1)/2)
    for i in range(len(label_list)):
        if len(label_list[i])>label_lengthest_size:
            label_lengthest_size = len(label_list[i])
    label_lengthest_size = int((label_lengthest_size+1)/2)

    batch_data_array = np.zeros(shape=(batch_size,data_lengthest_size,21),dtype=np.float32)
    batch_label_array = np.zeros(shape=(batch_size,label_lengthest_size,21),dtype=np.float32)
    for i in range(len(data_list)):
        index = data_list[i].split(",")
        for j in range(len(index)):
            batch_data_array[i][j][int(index[j])] = 1.0
    for i in range(len(label_list)):
        index = label_list[i].split(",")
        for j in range(len(index)):
            batch_label_array[i][j][int(index[j])] = 1.0
    return batch_data_array,batch_label_array

def get_rep_vec(file_path,vector_length):
    f = open(file_path)
    id_list = []
    seq_list = []
    while True:
        id = f.readline()
        if id == "":
            f.close()
            break
        seq = f.readline()
        id_list.append(id)
        seq_list.append(seq)
    rep_vec_array = np.zeros(shape=(len(seq_list),vector_length),dtype=np.float32)
    for i in range(len(seq_list)):
        index = seq_list[i].split(",")
        for j in range(len(index)):
            rep_vec_array[i][j] = float(index[j])
    return rep_vec_array,id_list

def get_Euclidean_Distance(result,rep_vec_array):
    result_list = []
    for i in range(rep_vec_array.shape[0]):
        result_list.append(np.sum(np.square(result - rep_vec_array[i])))
    return np.array(result_list, dtype=np.int32)
