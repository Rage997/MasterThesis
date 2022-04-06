import pandas as pd

'''
This package contains function for handling the dataset, filtering
with revelant data and reporting information about it
'''

def read_data(filename: str)->pd.DataFrame:
    df_orig = pd.read_excel(filename, sheet_name=None)
    # There are 3 sheets of the excel file. Get the correct one
    df = df_orig['GlobalAlienSpeciesFirstRecordDa']
    return df

def filter_data(df: pd.DataFrame)->pd.DataFrame:
    df = df[df['LifeForm'] == 'Viruses']
    df = df[df['FirstRecord'] > 1950]
    return df

def species_info(df: pd.DataFrame):
    family = df['LifeForm'].unique()
    region = df['Region'].unique()

    print(f'There are {len(family)} families and {len(region)} regions')
    species = df['TaxonName'].unique()

    print(f'There are {len(species)} species.')
    for p in family:
        df_s = df[df['LifeForm'] == p]
        animals = df_s['TaxonName'].unique()
        print(f'{p}[{len(animals)}]', end=', ')
