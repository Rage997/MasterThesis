# Load matrix into CUDA and see how much space it's using
import cupy as cp
mempool = cp.get_default_memory_pool()
print(mempool.used_bytes())              # 0
print(mempool.total_bytes())             # 0
M_gpu = cp.array(M)
est_bytes = np.prod(M.shape)*8
print(f'My estimate {est_bytes}') # num bytes
print(M.nbytes)                          # ~ num bytes
print(mempool.used_bytes())              # ~ num bytes
print(mempool.total_bytes())             # ~ num bytes

print(f'In MB {est_bytes/(10**6)}')
