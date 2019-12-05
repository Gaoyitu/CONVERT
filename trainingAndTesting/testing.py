import tensorflow as tf
import utils
import numpy as np
import copy
import time


residue_to_int = {'a':0,'r':1,'n':2,'d':3,'c':4,'q':5,'e':6,'g':7,'h':8,'i':9,'l':10,'k':11,'m':12,
                  'f':13,'p':14,'s':15,'t':16,'w':17,'y':18,'v':19,'<PAD>':20,'<UNK>':21,'<GO>':22,'<EOS>':23}
int_to_residue  = {idx: residue for residue, idx in residue_to_int.items()}

def get_data(test_file_path):
    data_list,target_list = utils.get_test_data(test_file_path)

    data_int = [[residue_to_int.get(residue,residue_to_int["<UNK>"])
                   for residue in line] for line in data_list]
    return data_int,target_list

batch_size = 1
gru_size = 300
embedding_size = 15
learning_rate = 1e-4
keep_prob = tf.placeholder(tf.float32,name="keep_prob")
inputs = tf.placeholder(tf.int32,[None,None],name="inputs")

def get_cell_fw():
    gru_cell_1 = tf.nn.rnn_cell.GRUCell(num_units=gru_size,
                                        kernel_initializer=tf.random_uniform_initializer(-0.1,
                                                                                         0.1,
                                                                                         seed=2))
    gru_cell_1_prob = tf.nn.rnn_cell.DropoutWrapper(cell=gru_cell_1, input_keep_prob=keep_prob, output_keep_prob=1.0)
    gru_cell_2 = tf.nn.rnn_cell.GRUCell(num_units=gru_size,
                                        kernel_initializer=tf.random_uniform_initializer(-0.1,
                                                                                         0.1,
                                                                                         seed=2))
    gru_cell_2_prob = tf.nn.rnn_cell.DropoutWrapper(cell=gru_cell_2, input_keep_prob=keep_prob, output_keep_prob=1.0)
    m_cell = tf.nn.rnn_cell.MultiRNNCell([gru_cell_1_prob, gru_cell_2_prob])
    return m_cell

def get_cell_bw():
    gru_cell_1 = tf.nn.rnn_cell.GRUCell(num_units=gru_size,
                                        kernel_initializer=tf.random_uniform_initializer(-0.1,
                                                                                         0.1,
                                                                                         seed=2))
    gru_cell_1_prob = tf.nn.rnn_cell.DropoutWrapper(cell=gru_cell_1, input_keep_prob=keep_prob, output_keep_prob=1.0)
    gru_cell_2 = tf.nn.rnn_cell.GRUCell(num_units=gru_size,
                                        kernel_initializer=tf.random_uniform_initializer(-0.1,
                                                                                         0.1,
                                                                                         seed=2))
    gru_cell_2_prob = tf.nn.rnn_cell.DropoutWrapper(cell=gru_cell_2, input_keep_prob=keep_prob, output_keep_prob=1.0)
    m_cell = tf.nn.rnn_cell.MultiRNNCell([gru_cell_1_prob, gru_cell_2_prob])
    return m_cell
def bi_encoder():
    encoder_embed_input = tf.contrib.layers.embed_sequence(inputs,len(residue_to_int),embedding_size)
    bi_encoder_output, bi_encoder_state = tf.nn.bidirectional_dynamic_rnn(
        cell_fw=get_cell_fw(),
        cell_bw=get_cell_bw(),
        inputs=encoder_embed_input,
        dtype=encoder_embed_input.dtype)
    encoder_output = tf.concat(bi_encoder_output, -1)
    encoder_state = []
    for layer_id in range(2):
        encoder_state.append(bi_encoder_state[0][layer_id])
        encoder_state.append(bi_encoder_state[1][layer_id])
    encoder_state = tuple(encoder_state)
    return encoder_output, encoder_state

def get_source(sources):
    for i in range(len(sources)):
        sources_i = sources[i]
        sources_i = np.array([sources_i],dtype=np.int32)
        yield sources_i

def start_test(input_file,output_file,rep_vec_file,net,vector_length,max_top_n):
    data_int,target_list = get_data(input_file)

    _, result = bi_encoder()
    output_f = open(output_file,"w")
    with tf.Session() as sess:
        rep_vec_array,id_list = utils.get_rep_vec(rep_vec_file,vector_length)
        saver = tf.train.Saver()
        sess.run(tf.global_variables_initializer())
        saver.restore(sess,net)
        for i, (data) in enumerate(get_source(sources=data_int)):
            time_start = time.time()
            rep_id_list = copy.deepcopy(id_list)
            output_f.write(target_list[i])
            result_value = sess.run(result, feed_dict={inputs: data,keep_prob:1.0})
            result_value = np.array(result_value).reshape(vector_length,)
            score_list = utils.get_Euclidean_Distance(result_value,rep_vec_array)
            for j in range(max_top_n):
                min_index = np.where(score_list == score_list.min())[0][0]
                if j == max_top_n - 1:
                    output_f.write(rep_id_list[min_index])
                else:
                    output_f.write(rep_id_list[min_index].strip("\n") + ',')

                score_list = np.delete(score_list, [min_index])
                del rep_id_list[min_index]
            time_end = time.time()
            print("iterater:",i,"time:",time_end-time_start)

if __name__ =="__main__":
    start_test("data/test/SCOP95.txt",
               "data/result/SCOP.txt",
               "data/vector/SCOP95.txt",
               "data/net/SCOP/SCOP.ckpt",
               1200,
               100)



