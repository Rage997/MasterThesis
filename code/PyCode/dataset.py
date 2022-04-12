from time import time
import pandas as pd
import numpy as np

'''
This package contains function for handling the dataset, filtering
with revelant data, reporting information about it and converting it
to a matrix to be lately processed with a REM model
'''
class InvasiveSpecies:
    def __init__(self, filename: str) -> None:
       
        try:
            self.df = self._read_data(filename)
        except FileNotFoundError as e:
            raise FileNotFoundError('File not found. Check data path in config.py')
        
        self.filter_data()

        self.species = self.df['TaxonName'].unique()
        self.region = self.df['Region'].unique()
        self.Ns = len(self.species)
        self.Nr = len(self.region)

    def _read_data(self, filename: str)->pd.DataFrame:
        df_orig = pd.read_excel(filename, sheet_name=None)
        # There are 3 sheets of the excel file. Get the correct one
        df = df_orig['GlobalAlienSpeciesFirstRecordDa']
        return df

    def get_time_interval(self):
        t_min = self.df['FirstRecord'].min()
        t_max = self.df['FirstRecord'].max()
        return (t_min, t_max)

    def filter_data(self)->None:
        # self.df = self.df[self.df['LifeForm'] == 'Viruses']
        self.df = self.df[self.df['FirstRecord'] > 1950]
        # return df

    def remove_irrelevant(self, species_tol: int, region_tol: int):
        '''Removes irrelevant nodes'''

        df = self.df
        spec_iter = np.zeros(len(self.species), dtype=int)
        reg_iter = np.zeros(len(self.region), dtype=int)

        for idx, s in enumerate(self.species):
            df_s = df[df['TaxonName'] == s]
            # n_iter = len(df_s)
            spec_iter[idx] = len(df_s)
            # print(f'The species {s} has {spec_iter[idx]} iteractions')

        for idx, s in enumerate(self.region):
            df_s = df[df['Region'] == s]
            # n_iter = len(df_s)
            reg_iter[idx] = len(df_s)
            # print(f'The species {s} has {spec_iter[idx]} iteractions')

        # select the species that have less than 3 iteractions
        remove_s_idx = np.where(spec_iter < species_tol)
        # and regions that are invaded by less than 10 species
        remove_r_idx = np.where(reg_iter < region_tol)

        # Filter the dataset: remove all rows belonging to species with less than 3 iteractions
        for s in self.species[remove_s_idx]:
            # print(s)  
            df = df.drop(df[df['TaxonName'] == s].index)

        # Filter the dataset: remove all rows belonging to region with less than 10 invasions
        for s in self.region[remove_r_idx]:
            # print(s)  
            df = df.drop(df[df['Region'] == s].index)

        species_removed = len(self.species) - len(df['TaxonName'].unique())
        print(f'Removed {species_removed} species')
        region_removed = len(self.region) - len(df['Region'].unique())
        print(f'Removed {region_removed} region')
        
        # update class variables
        self.region = df['Region'].unique()
        self.species = df['TaxonName'].unique()
        self.df = df

    def log_species_info(df: pd.DataFrame)->None:
        family = df['LifeForm'].unique()
        region = df['Region'].unique()

        print(f'There are {len(family)} families and {len(region)} regions')
        species = df['TaxonName'].unique()

        print(f'There are {len(species)} species.')
        for p in family:
            df_s = df[df['LifeForm'] == p]
            animals = df_s['TaxonName'].unique()
            print(f'{p}[{len(animals)}]', end=', ')

    def build_matrix(self, time_resolution:int=5)->np.array:

        t_min, t_max = self.get_time_interval()
        self.time = np.arange(t_min, t_max, time_resolution)

        # Non square matrix M \in R^{T, Ns, Nr}
        M = np.zeros((len(self.time), self.Ns, self.Nr))
        # or...build full matrix
        # p = self.Ns + self.Nr
        # M = np.zeros((len(self.time), p, p)) # if full
    
        for i, t in enumerate(self.time):
            df_now = self.df[(self.df['FirstRecord'] >= t) & (self.df['FirstRecord'] < t+(time_resolution-1))]
            for _, row in df_now.iterrows():
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