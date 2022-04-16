from re import L
import pandas as pd
import numpy as np
# import scipy
# from scipy.spatial import distance_matrix
# see: https://sparrow.dev/pairwise-distance-in-numpy/
# may be useful if we want to run some code on gpu
from scipy.spatial.distance import cdist

# TODO remove this import and refactor the class to reduce coupling
# from config import d
# TODO I honestly don't think that I need a class for this code.


np.random.seed(42)
class EKR:
    def __init__(self, Y_kf: np.array) -> None:
        self.Y_kf = Y_kf
        self.Nt, self.Ns, self.Nr  = Y_kf.shape
        print(f"Matrix dim is {Y_kf.shape}")
        self.p_y = self.Ns * self.Nr
        self.dim = 2 # latent space dimension
        
        alpha = 1
        self.gam_sum_mat = np.full((self.Nt, self.p_y), np.exp(alpha))

    def h(self, X: np.array, gam: np.array):
        # Compute pairwise distance
        X = X.reshape((-1, self.dim))
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
        H = np.zeros((self.p_y, self.p_y*self.self.dim))
        ii = 0
        for s in range(self.Ns):
            for r in range(self.Ns+1, self.p_y):
                # TODO not sure if it's self.dim-1 or just self.dim, check it
                sub_idx1 = self.dim*s - np.arange(self.dim-1, -1, -1) # index species
                sub_idx2 = self.dim*r - np.arange(self.dim-1, -1, -1) # index regions
                query = np.hstack((sub_idx1, sub_idx2))
                H[ii, query] = self.d_h(x_prior[sub_idx1], x_prior[sub_idx2], gam[ii])
                ii += 1
            
    def run(self, n_iter: int, past_censoring=False):
        
        pd = self.p_y*self.dim
        x_prior = x_post = np.zeros((self.Nt, pd))
        p_prior = p_post = np.zeros((self.Nt, pd, pd))
        Q = np.diag(np.full(pd, 0.001))

        P0_smooth = P_0 = Q # Starting points
        x0_smooth = x0 = np.random.randn(pd)

        self.build_H(x_prior, 0) #TODO

        for i in range(n_iter):
            for t in range(self.Nt):
                if t == 0: # init, see pg 9 igor paper
                    x_prior[t, :] = x0
                    p_prior[t] = P_0 + Q
                else: # a-prior distribution
                    x_prior[t] = x_post[t-1]
                    p_prior[t] = p_post[t-1] + Q
                    
                # see pg 9 igor paper
                H = self.build_H(x_prior[t], self.gam_sum_mat[t]) # TODO
                # you can do taylor: mu(x) = mu(x0) + d_h(x - x0)
                mu = self.h(x_prior[t], self.gam_sum_mat[t])
                if (t>1 and past_censoring):
                    mu = mu * 0 # TODO censor past
                
                vv = mu
                R = np.diag(vv)
                R_inv = 1/R
                R_inv[R_inv == np.inf] = 0

                # Sherman-Morrison-Woodbury identity (page 31) (kinda)
                # K = np.linalg.solve(p_prior[t])
