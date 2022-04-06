import pandas as pd
import numpy as np

class EKR:
    def __init__(self, dataset: pd.DataFrame) -> None:
        # TODO there's a high coupling with the dataset structure here
        # If I want to make it more general I should replace data with 
        # sender and receiver or something like it.

        self.df = dataset
        self.species = self.df['TaxonName'].unique()
        self.region = self.df['Region'].unique()
        self.Ns = len(self.species)
        self.Nr = len(self.region)

        t_min = self.df['FirstRecord'].min()
        t_max = self.df['FirstRecord'].max()
        self.time = np.arange(t_min, t_max, 2)

    # TODO get methods for:
    # Get specie name by id and viceversa
    # def get_spec_idx(name: str)->int:
        

    def build_matrix(self)->np.array:
        
        # species = df['TaxonName'].unique()
        # region = df['Region'].unique()

        # # Buld matrix
        # n_s = len(species)
        # n_r = len(region)

        # t_min = df['FirstRecord'].min()
        # t_max = df['FirstRecord'].max()
        # time = np.arange(t_min, t_max, 2)

        p = self.Ns + self.Nr
        M = np.zeros((len(self.time), self.Ns, self.Nr))
        # M = np.zeros((len(time), p, p)) # if full

        for i, t in enumerate(self.time):
            df_now = self.df[(self.df['FirstRecord'] >= t) & (self.df['FirstRecord'] < t+1)]
            for index, row in df_now.iterrows():
                s = row['TaxonName']
                r = row['Region']
                # print(f'Species {s} invaded region {r} at time {t}')

                s_idx = np.where(self.species == s)
                r_idx = np.where(self.region == r)
                M[i, s_idx, r_idx] = 1
        print(M.shape)
        return M

    def export_to_R(self):
        # Export data and then import it into R
        # print(n_s, n_r)
        M = M.reshape(len(self.time), self.Ns*self.Nr)
        # M = M.reshape(len(time), (p)**2) # if full matix
        # print(f's = {n_s}, r = {n_r}')
        # print(M.shape)
        np.save('matrix_full.npy', M)