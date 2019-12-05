import tensorflow as tf
import utils
import numpy as np

residue_to_int = {'a': 0, 'r': 1, 'n': 2, 'd': 3, 'c': 4, 'q': 5, 'e': 6, 'g': 7, 'h': 8, 'i': 9, 'l': 10, 'k': 11,
                      'm': 12,
                      'f': 13, 'p': 14, 's': 15, 't': 16, 'w': 17, 'y': 18, 'v': 19, '<PAD>': 20, '<UNK>': 21, '<GO>': 22,
                      '<EOS>': 23}
int_to_residue = {idx: residue for residue, idx in residue_to_int.items()}

def get_data(test_file_path):
    data_list, target_list = utils.get_test_data(test_file_path)

    data_int = [[residue_to_int.get(residue, residue_to_int["<UNK>"])
                 for residue in line] for line in data_list]
    return data_int,target_list

batch_size = 1
gru_size = 300
embedding_size = 15
learning_rate = 1e-4
keep_prob = tf.placeholder(tf.float32, name="keep_prob")
inputs = tf.placeholder(tf.int32, [None, None], name="inputs")


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
    encoder_embed_input = tf.contrib.layers.embed_sequence(inputs, len(residue_to_int), embedding_size)
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
        sources_i = np.array([sources_i], dtype=np.int32)
        yield sources_i

def start_generate(input_file,output_file,net_file):

    _,encoder_state = bi_encoder()
    data_int,id_list = get_data(input_file)
    with tf.Session() as sess:
        saver = tf.train.Saver()
        sess.run(tf.global_variables_initializer())
        saver.restore(sess,net_file)
        write = open(output_file, "w")

        for i, (data) in enumerate(get_source(sources=data_int)):
            a_value = sess.run(encoder_state, feed_dict={inputs: data, keep_prob: 1.0})
            a_array = np.array(a_value)
            write.write(id_list[i])
            for index in range(4):
                for axis_3 in range(300):
                    if index == 3:
                        if axis_3 == 299:
                            write.write(str(a_array[index][0][axis_3]) + "\n")
                        else:
                            write.write(str(a_array[index][0][axis_3]) + ",")
                    else:
                        write.write(str(a_array[index][0][axis_3]) + ",")
            print(i)

if __name__ =="__main__":
    start_generate("finalData/train/SCOP.txt",
                   "finalData/vector/SCOP.txt",
                   "finalData/net/SCOP/SCOP.ckpt")



