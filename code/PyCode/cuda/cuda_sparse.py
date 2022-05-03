import cupy as cp
from cupyx.scipy.sparse import csr_matrix as csr_gpu
from scipy.sparse import random

sp = random(22, n_s*n_r, density=0.25)
sp_gpu = csr_gpu(sp)
x = np.random.rand(n_s*n_r)

res = sp.dot(x)

print(sp_gpu.shape, x.shape)
x_gpu = cp.array(x) #moving x to the gpu
print(sp_gpu.shape[1] != len(x_gpu))

res_gpu = sp_gpu.dot(x_gpu)

if np.allclose(res, cp.asnumpy(res_gpu)):
    print("They are the same!")
else:
    print("They are different!")
