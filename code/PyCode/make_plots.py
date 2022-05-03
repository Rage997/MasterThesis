from dataset import InvasiveSpecies
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def histogram(data: InvasiveSpecies):
    count_animals = []
    for p in data.families:
        df_s = data.df[data.df['LifeForm'] == p]
        animals = df_s['TaxonName'].unique()
        animals = len(animals)
        count_animals.append(animals)
        print(f'{p}[{animals}]', end=', ')

    plt.bar(range(len(count_animals)), count_animals, align='center')
    firstletter_family = [f[:3]+'.' for f in data.families]
    plt.xticks(range(len(data.families)), firstletter_family, size='small')
    plt.savefig('../../latex/figures/histogram_taxfam.png', dpi=1200)
    plt.show()


def total_average_invasion_per_year(data: InvasiveSpecies):

    avg_inv = np.zeros(len(data.time))
    for idx in range(len(data.time)):
        avg_inv[idx] = np.sum(data.M[idx])

    plt.plot(data.time, avg_inv)
    plt.savefig('../../latex/figures/total_avg_invasion.png', dpi=1200)

    plt.show()

def species_region_invasion(data: InvasiveSpecies):
    '''Plots the amount of invasion for each species'''

    spec_inv = []
    # time_res = data.time_resolution
    # for i, t in enumerate(data.time):
    #     df_now = data.df[(data.df['FirstRecord'] >= t) & (data.df['FirstRecord'] < t+(time_res-1))]
    for s in data.species:
        df_s = data.df[data.df['TaxonName'] == s]
        invasion = df_s['Region'] # get all the invasion of species s
        invasion = len(invasion)
        spec_inv.append(invasion)
    
    # Display the probability density function (pdf)
    spec_inv.sort()
    hmean = np.mean(spec_inv)
    hstd = np.std(spec_inv)
    pdf = stats.norm.pdf(spec_inv, hmean, hstd)
    print(f'mean: {hmean}, std{hstd}')
    plt.plot(spec_inv, pdf) 
    plt.show()



