import numpy as np

def distance_matrix(x: np.array, y: np.array) -> np.array:
    # TODO does this work on 3D 
    np.linalg.norm(x[:, None, :] - y[None, :, :], axis=-1)

