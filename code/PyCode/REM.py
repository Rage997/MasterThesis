import pandas as pd
import numpy as np
# import scipy
# from scipy.spatial import distance_matrix
# see: https://sparrow.dev/pairwise-distance-in-numpy/
# may be useful if we want to run some code on gpu
from scipy.spatial.distance import cdist

# TODO remove this import and refactor the class to reduce coupling
from config import d
class EKR:
    def __init__(self, Y_kf: np.array) -> None:
       self.Y_kf = Y_kf

    def h(self, X: np.array, gam: np.array):
        # Compute pairwise distance
        X = X.reshape((-1, d))
        D = cdist(X[:, 0], X[:, 1]) #TODO should I square?
        # Distances are symmetric -> get lower triangular
        # and flaten it to vector
        D = D[np.tril(D) != 0].flatten()
        return gam * np.exp(D)

    def d_h(self, xi, xj, gam_ij):  
        dist = np.linalg.norm(xi - xj)
        return 2*gam_ij * np.exp(-dist) * np.concatenate((xj - xi, xi - xj), axis=0)

    def build_H(self, x_prior: np.array, gam: np.array) -> np.array:
        # Populate H_ij = [derivative function(i,j)] -> sparse matrix
        H = np.zeros((self.p_y, self.p_y*self.self.d))
        ii = 0
        for s in range(self.Ns):
            for r in range(self.Ns+1, self.p_y):
                # TODO not sure if it's self.d-1 or just self.d, check it
                sub_idx1 = self.d*s - np.arange(self.d-1, -1, -1) # index species
                sub_idx2 = self.d*r - np.arange(self.d-1, -1, -1) # index regions
                query = np.hstack((sub_idx1, sub_idx2))
                H[ii, query] = self.d_h(x_prior[sub_idx1], x_prior[sub_idx2], gam[ii])
                ii += 1
            

    def run(self):
        x_prior = x_post = np.zeros((n, p*d))
        p_prior = p_post = np.zeros((p*d, p*d, n))
        Q = np.diag(np.full(p*d, 0.001))

        P0_smooth = P_0 = Q # Starting points
        x0_smooth = x0 = np.random.randn(p*d)

        self.build_H(x_prior, 0) #TODO

