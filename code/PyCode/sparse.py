import tensorflow as tf
import scipy.sparse as ss
import numpy as np


def scipy_run(sess, result, X, tf_indices, tf_id_values, tf_weight_values,
              tf_dense_shape):
    row_nnz = np.diff(X.indptr)
    indices = np.asarray([[row_i, col_i]
                          for row_i, nnz in enumerate(row_nnz)
                          for col_i in range(nnz)], dtype=np.int64)
    ids = X.indices.astype(np.int64)
    weights = X.data
    tf_result = sess.run(result, {tf_indices: indices,
                                  tf_id_values: ids,
                                  tf_weight_values: weights,
                                  tf_dense_shape: X.shape})
    return tf_result

def main():
    X = ss.csr_matrix([[1, 0, 0],
                       [1, 0, 1],
                       [0, 0, 1],
                       [0, 2, 0]], dtype=np.float32)
    W = np.asarray([[.1, .2, .3],
                    [.2, .3, .1],
                    [.3, .2, .2]], dtype=np.float32)
    b = np.asarray([[.1, .2, .3]], dtype=np.float32)

    # scipy version
    direct_result = X @ W + b

    print(direct_result)

    # tensorflow version

    def convert_sparse_matrix_to_sparse_tensor(X):
        coo = X.tocoo()
        indices = np.mat([coo.row, coo.col]).transpose()
        return tf.SparseTensor(indices, coo.data, coo.shape)


    # tf_indices = tf.keras.Input(dtype = tf.int64, shape=(None, 2))
    # tf_id_values = tf.keras.Input(dtype = tf.int64)
    # tf_weight_values = tf.keras.Input(dtype = tf.float32)
    # tf_dense_shape = tf.keras.Input(tf.int64, [2])
    # sp_ids = tf.SparseTensor(tf_indices, tf_id_values, tf_dense_shape)
    # sp_weights = tf.SparseTensor(tf_indices, tf_weight_values, tf_dense_shape)

    tf_X = convert_sparse_matrix_to_sparse_tensor(X)
    tf_W = tf.Variable(W, tf.float32)
    tf_b = tf.Variable(b, tf.float32)
    
    tf_result = X @ W + b

    print(tf_result)


    if (direct_result == tf_result).all():
        print("They are the same!")
    else:
        print("They are different!")


if __name__ == '__main__':
    main()