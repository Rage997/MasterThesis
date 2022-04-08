import pandas as pd
import numpy as np

'''
This package contains function for handling the dataset, filtering
with revelant data, reporting information about it and converting it
to a matrix to be lately processed with a REM model
'''

class InvasiveSpecies:
    def __init__(self, filename: str) -> None:
        # TODO there's a high coupling with the dataset structure here
        # If I want to make it more general I should replace data with 
        # sender and receiver or something like it.

        try:
            self.df = self._read_data(filename)
        except FileNotFoundError as e:
            raise FileNotFoundError('File not found. Check data path in config.py')
        
        self.filter_data()

        self.species = self.df['TaxonName'].unique()
        self.region = self.df['Region'].unique()
        self.Ns = len(self.species)
        self.Nr = len(self.region)

        t_min = self.df['FirstRecord'].min()
        t_max = self.df['FirstRecord'].max()
        self.time = np.arange(t_min, t_max, 2)

    def _read_data(self, filename: str)->pd.DataFrame:
        df_orig = pd.read_excel(filename, sheet_name=None)
        # There are 3 sheets of the excel file. Get the correct one
        df = df_orig['GlobalAlienSpeciesFirstRecordDa']
        return df

    def filter_data(self)->None:
        self.df = self.df[self.df['LifeForm'] == 'Viruses']
        self.df = self.df[self.df['FirstRecord'] > 1950]
        # return df

    def species_info(df: pd.DataFrame)->None:
        family = df['LifeForm'].unique()
        region = df['Region'].unique()

        print(f'There are {len(family)} families and {len(region)} regions')
        species = df['TaxonName'].unique()

        print(f'There are {len(species)} species.')
        for p in family:
            df_s = df[df['LifeForm'] == p]
            animals = df_s['TaxonName'].unique()
            print(f'{p}[{len(animals)}]', end=', ')

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
        # print(M.shape)
        self.M = M.copy() # store it into class
        return self.M

    def export_to_R(self):
        # Export data and then import it into R
        # print(n_s, n_r)
        M = self.M.reshape(len(self.time), self.Ns*self.Nr)
        # M = M.reshape(len(time), (p)**2) # if full matix
        print(f's = {self.Ns}, r = {self.Nr}')
        print(M.shape)
        np.save('../data/matrix.npy', M)