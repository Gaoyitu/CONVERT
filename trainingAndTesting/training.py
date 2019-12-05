import tensorflow as tf
from tensorflow.python.layers.core import Dense
import utils
import numpy as np

train_file_path = "data/train/SCOP.txt"
#train_file_path = "data/train/SCOPe.txt"
#train_file_path = "data/test.txt"
residue_to_int = {'a':0,'r':1,'n':2,'d':3,'c':4,'q':5,'e':6,'g':7,'h':8,'i':9,'l':10,'k':11,'m':12,
                  'f':13,'p':14,'s':15,'t':16,'w':17,'y':18,'v':19,'<PAD>':20,'<UNK>':21,'<GO>':22,'<EOS>':23}
int_to_residue  = {idx: residue for residue, idx in residue_to_int.items()}

data_list,target_list = utils.get_shuffle_data(train_file_path)
data_int = [[residue_to_int.get(residue,residue_to_int["<UNK>"])
               for residue in line] for line in data_list]
target_int = [[residue_to_int.get(residue,residue_to_int["<UNK>"])
              for residue in line]+[residue_to_int["<EOS>"]] for line in target_list]

epoch_size = 200
batch_size = 32
gru_size = 300
embedding_size = 15
learning_rate = 1e-4
keep_prob = tf.placeholder(tf.float32,name="keep_prob")
inputs = tf.placeholder(tf.int32,[None,None],name="inputs")
targets = tf.placeholder(tf.int32,[None,None],name="targets")
target_sequence_length = tf.placeholder(tf.int32,(None,),name="target_sequence_length")
max_target_sequence_length = tf.reduce_max(target_sequence_length,name="max_target_len")

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

def attention_decoder_cell(encoder_output):
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
    gru_cell_3 = tf.nn.rnn_cell.GRUCell(num_units=gru_size,
                                        kernel_initializer=tf.random_uniform_initializer(-0.1,
                                                                                         0.1,
                                                                                         seed=2))
    gru_cell_3_prob = tf.nn.rnn_cell.DropoutWrapper(cell=gru_cell_3, input_keep_prob=keep_prob, output_keep_prob=1.0)
    gru_cell_4 = tf.nn.rnn_cell.GRUCell(num_units=gru_size,
                                        kernel_initializer=tf.random_uniform_initializer(-0.1,
                                                                                         0.1,
                                                                                         seed=2))
    gru_cell_4_prob = tf.nn.rnn_cell.DropoutWrapper(cell=gru_cell_4, input_keep_prob=keep_prob, output_keep_prob=1.0)
    m_cell = tf.nn.rnn_cell.MultiRNNCell([gru_cell_1_prob, gru_cell_2_prob,gru_cell_3_prob,gru_cell_4_prob])
    attention_mechanim = tf.contrib.seq2seq.BahdanauAttention(gru_size,
            encoder_output, normalize = True)
    cell = tf.contrib.seq2seq.AttentionWrapper(m_cell, attention_mechanim,
            attention_layer_size=gru_size)
    return cell
def decoder_projection():
    output_layer = Dense(units=len(residue_to_int),
                         kernel_initializer=tf.truncated_normal_initializer(mean=0.0, stddev=0.1))
    return output_layer
def process_decoder_input():
    ending = tf.strided_slice(targets,[0,0],[batch_size,-1],[1,1])
    decoder_input = tf.concat([tf.fill([batch_size,1],residue_to_int["<GO>"]),ending],1)
    return decoder_input
def train_decoder(encoder_output,encoder_state):
    target_vocab_size = len(residue_to_int)
    decoder_embeddings = tf.Variable(tf.random_uniform([target_vocab_size,embedding_size]))
    decoder_input = process_decoder_input()
    decoder_embed_input = tf.nn.embedding_lookup(decoder_embeddings,decoder_input)
    projection_layer = decoder_projection()
    decoder_cell = attention_decoder_cell(encoder_output)

    init_state = decoder_cell.zero_state(batch_size, tf.float32).clone(cell_state=encoder_state)
    helper = tf.contrib.seq2seq.TrainingHelper(decoder_embed_input,
                                               target_sequence_length,
                                               time_major=False)
    decoder = tf.contrib.seq2seq.BasicDecoder(decoder_cell,
                                              helper,
                                              init_state,
                                              output_layer=projection_layer)
    outputs, _, _ = tf.contrib.seq2seq.dynamic_decode(decoder,
                                                      impute_finished  = True,
                                                      maximum_iterations=max_target_sequence_length)
    return outputs.rnn_output
def seq2seq_model():
    encoder_output, encoder_state = bi_encoder()
    decoder_output = train_decoder(encoder_output,encoder_state)
    return decoder_output
def pad_sentence_batch(sentence_batch):
    max_sentence = max([len(sentence) for sentence in sentence_batch])
    return [sentence+[residue_to_int["<PAD>"]]*(max_sentence - len(sentence)) for sentence in sentence_batch]
def get_batches(sources,targets):
    for batch_i in range(0,len(sources)//batch_size):
        start_i = batch_i*batch_size
        sources_batch = sources[start_i:start_i+batch_size]
        targets_batch = targets[start_i:start_i + batch_size]
        pad_sources_batch = np.array(pad_sentence_batch(sources_batch))
        pad_targets_batch = np.array(pad_sentence_batch(targets_batch))
        targets_lengths = []
        for target in targets_batch:
            targets_lengths.append(len(target))
        yield pad_sources_batch,pad_targets_batch,targets_lengths
def cost():
    training_logits = seq2seq_model()
    masks = tf.sequence_mask(target_sequence_length,max_target_sequence_length,dtype=tf.float32,name="mask")
    loss = tf.contrib.seq2seq.sequence_loss(training_logits,targets, masks)
    return loss


def start_train(net):
    loss = cost()
    train_step = tf.train.AdamOptimizer(learning_rate).minimize(loss)
    saver = tf.train.Saver()
    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        #saver.restore(sess,net)
        for epoch in range(epoch_size):
            for batch_i, (data_batch, target_batch,targets_lengths) in enumerate(
                    get_batches(data_int,
                                target_int)):
                sess.run(train_step, feed_dict={inputs: data_batch,
                                                 targets:target_batch,
                                                 target_sequence_length:targets_lengths,
                                                 keep_prob:0.7})
                if batch_i%2 == 0:
                    loss_value = sess.run(loss,feed_dict={inputs: data_batch,
                                                     targets:target_batch,
                                                     target_sequence_length:targets_lengths,
                                                     keep_prob:1.0})
                    print("epoch:",epoch+1,"batch:",batch_i+1,"loss:",loss_value)
        saver.save(sess,net)
if __name__ == "__main__":
    start_train("data/net/SCOP/SCOP.ckpt")
    #start_train("data/net/SCOPe/SCOPe.ckpt")


