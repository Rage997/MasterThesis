# simulate a clustering dynamics process
import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
import numpy as np
import pandas as pd
from tensorflow.python.client import device_lib
import matplotlib.pyplot as plt


# Hide GPU from visible devices
# tf.config.set_visible_devices([], 'GPU')

gpus = tf.config.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)


# -------------------- Load real data -------------------------
# Load data from memory
Y = np.load('../data/result/452_75/matrix.npy')
d = tf.constant(2, dtype =tf.int32)
n = tf.constant(Y.shape[0], dtype =tf.int32)
s = tf.constant(Y.shape[1], dtype =tf.int32) 
r = tf.constant(Y.shape[2], dtype =tf.int32)

Y = Y.reshape(n, s*r)
Y = tf.Variable(Y, dtype = tf.float32) #flatten it to array

n_nodes = tf.constant(s+r, dtype =tf.int32)
print(f'There are {n_nodes} nodes')
time_interval = tf.constant(range(1, n+1), dtype =tf.int32)
X_true = np.zeros((n, n_nodes, d))
alpha = tf.constant(3, dtype = tf.float32)
Z = tf.constant(0., dtype = tf.float32)
for i in range(n_nodes):
  for j in range(d):
    X_true[ :, i, j] = (tf.keras.activations.sigmoid(np.random.uniform(-2, 2, 1)*10*(time_interval/n-np.random.uniform(0.2, 0.8, 1))) -0.5)*np.random.uniform(1, 3, 1) + np.random.uniform(-0.5, 0.5, 1)

X_true = tf.constant(X_true, dtype = tf.float32)


# ------------------- Model declaration -------------

# def h_b(X_b, b, c):
#   species = tf.reshape(X_b[:, :(s*d)], (b, s, d))
#   region = tf.reshape(X_b[:, (s*d):], (b, r, d))
#   return tf.exp( c + tf.reshape(tf.matmul(species, region, transpose_b=True), (b, s*r)))

def h(x, c):
  species = tf.reshape(x[:(s*d)], (s, d))
  region = tf.reshape(x[(s*d):], (r, d))
  return tf.exp( c + tf.reshape(tf.matmul(species, region, transpose_b=True), [-1]))

def get_offset(x):
  species = tf.reshape(x[:(s*d)], (s, d))
  region = tf.reshape(x[(s*d):], (r, d))
  return tf.reshape(tf.matmul(species, region, transpose_b=True), [-1])

@tf.function
def get_h_H(x, c): 
  with tf.GradientTape(persistent=True) as tape:
    tape.watch(x)
    h_current = h(x, c)
  H = tape.jacobian(target = h_current, sources = x, parallel_iterations=1000 , experimental_use_pfor=True)
  return h_current, H


def get_censoring(Y):
  censoring = tf.Variable(np.ones((n, s*r)), dtype = tf.float32)
  for t in range(1,n):
    censoring[t, :].assign(censoring[t-1, :]-Y[t-1, :])
  return tf.keras.activations.relu(censoring)

class Kalman_model(tf.Module):
  def __init__(self, x_0, sigma2):
    #super().__init__(**kwargs)
    self.x = tf.Variable(x_0)
    self.V = tf.Variable(sigma2)
    #self.H = tf.Variable(tf.zeros((s*r, n_nodes*d)))

  def __call__(self, y, c, censoring):
    self.V.assign_add(sigma2)
    mu, H = get_h_H(self.x, c)
    # X_b =  tf.reshape(self.x, [b, n_nodes*d]) + tf.sqrt(sigma2) * tf.random.normal(shape = [b, n_nodes*d], dtype = tf.float32)
    # Y_b = h_b()
    R_inv = 1/mu
    if(censoring != None):
      censoring = tf.squeeze(censoring)
      H = tf.transpose(censoring * tf.transpose(H))
      # H = tf.matmul(tf.linalg.diag(censoring), H)
      R_inv = 1/mu * censoring
      mu = mu * censoring
    S_chol = tf.linalg.cholesky(tf.linalg.diag(1/self.V) + tf.matmul(H, tf.transpose(R_inv * tf.transpose(H)), transpose_a=True))
    K = tf.linalg.cholesky_solve(S_chol, R_inv * tf.transpose(H))
    #K = tf.matmul(tf.linalg.inv(tf.linalg.diag(1/self.V) + tf.matmul(H, tf.matmul(R_inv, H), transpose_a=True)), tf.matmul(H, R_inv, transpose_a=True))
    self.x.assign_add( tf.squeeze(tf.matmul(K, tf.reshape(y - mu, (-1,1)))))
    self.V.assign((tf.cast(1, dtype=tf.float32) - tf.reduce_sum(K*tf.transpose(H), axis=1)) *self.V)
    return 

  def smoother(self, x_prev, V_prev):
    V_prior = V_prev + sigma2
    B = V_prev/V_prior
    self.x.assign(x_prev + tf.squeeze(B * tf.squeeze(self.x - x_prev)))
    self.V.assign( V_prev + tf.square(B)*(self.V - V_prior))
    return 

#train_dataset_rev = tf.data.Dataset.from_tensor_slices(tf.reverse(Y, axis=[0]))
#train_dataset_rev = train_dataset.batch(batch_size=1)

@tf.function
def train_kalman(model):
  # tf.print("forward filtering")
  for t, (y_batch, censor_batch) in enumerate(train_dataset):
    model(y_batch, c[t, :], censor_batch)
    X_kalman[t+1, :].assign(model.x)
    V_kalman[t+1, :].assign(model.V)
    # tf.print(t)


@tf.function
def train_smoother(model):
  # tf.print("backward smoothing")
  for t, (y_batch, censor_batch) in enumerate(train_dataset):
    t = n-tf.cast(t, dtype=tf.int32)
    if(t!=0):
      offset[t-1, :].assign(get_offset(model.x))
    model.smoother(X_kalman[t-1, :], V_kalman[t-1, :])
    X_kalman[t-1, :].assign(model.x)
    V_kalman[t-1, :].assign(model.V)
    # tf.print(t)

@tf.function(autograph=False)
def fit_glm(z, y, offset, censor_data):
  if(censor_data != None):
    wh = tf.squeeze(tf.where(tf.reshape(censor_data, [-1])==1))
    z = tf.gather(z, wh, axis=0)
    y = tf.gather(y, wh, axis=0)
    offset = tf.gather(offset, wh, axis=0)
  model_coefficients, linear_response, is_converged, num_iter = tfp.glm.fit(model_matrix=z, response=y, model=tfp.glm.Poisson(), offset = offset )
  log_likelihood = tfp.glm.Poisson().log_prob(y, linear_response)
  return (model_coefficients, log_likelihood)


# censor_data = None
censor_data = get_censoring(Y)

train_dataset = tf.data.Dataset.from_tensor_slices((Y, censor_data))
train_dataset = train_dataset.batch(batch_size=1)

# for t, (y_batch, censor_batch) in enumerate(train_dataset):
#   print(censor_batch)
#   print(y_batch)

model_coefficients = alpha
c = tf.Variable(tf.fill((n, s*r), model_coefficients), trainable=False)
X_kalman = tf.Variable(tf.random.uniform(shape = (n+1,  n_nodes * d) , minval= -1, maxval= 1, dtype=tf.float32))
V_kalman = tf.Variable(tf.random.uniform(shape = (n+1,  n_nodes * d) , minval= -1, maxval= 1, dtype=tf.float32))
offset = tf.Variable(tf.random.uniform(shape = (n, s*r ) , minval= -1, maxval= 1, dtype=tf.float32))
Z = tf.ones((n, s*r), dtype=tf.float32)

#x_0 = tf.random.uniform(shape = [n_nodes * d] , minval= -1, maxval= 1, dtype=tf.float32)
x_0 = tf.reshape(X_true[0, :, :], -1)
sigma2 = tf.constant( [0.001]*(n_nodes*d).numpy(), dtype=tf.float32)

X_kalman[0,:].assign(x_0)
V_kalman[0,:].assign(sigma2)
 
# -------------------- Run model -------------------------

model = Kalman_model(x_0, sigma2)
#tf.config.run_functions_eagerly(False)
for iter in range(50):
  train_kalman(model)
  train_smoother(model)
  # model_coefficients, log_likelihood = fit_glm(tf.reshape(Z, (-1,1)), tf.reshape(Y, [-1]), tf.reshape(offset, [-1]), censor_data)
  model_coefficients = tf.math.log(tf.reduce_sum(Y)/tf.reduce_sum(tf.math.exp(offset)*censor_data))
  c.assign(tf.fill((n, s*r), model_coefficients))
  # tf.print("log likelihood = ", tf.reduce_sum(log_likelihood))
  tf.print(f"Interation number {i}")
  tf.print("alpha = ", model_coefficients)


# --------------------- Run plots ---------------------
Y = tf.reshape(Y, (n, s*r))
Y_est = tf.Variable(np.zeros(shape = (n, s*r)), dtype=tf.float32)
for t in range( n):
   Y_est[t, :].assign(h(X_kalman[t+1, :], c[t, :]))


n_plots = 10
i = tf.random.uniform(shape = [n_plots], maxval = s, dtype = tf.int32)
j = tf.random.uniform(shape = [n_plots], maxval = r, dtype = tf.int32)
fig, axs = plt.subplots(n_plots, n_plots, figsize=(25,25))
for i in range(n_plots):
  for j in range(n_plots):
    axs[ i, j].plot( Y[:, i + j*r ])
    axs[ i, j].plot( Y_est[:, i + j*r ])

path = f'../data/result/{s}_{r}/'
plt.savefig(path+"model_result.png", dpi=75)

# Save results
tf.print("alpha = ", model_coefficients)
np.save(path+'resultX', X_kalman)
np.save(path+'resultY', Y_est)
