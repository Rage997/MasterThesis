from dataset import InvasiveSpecies
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

plot_dir = '../../latex/figures/'

def histogram(data: InvasiveSpecies):
    count_animals = []
    for p in data.families:
        df_s = data.df[data.df['LifeForm'] == p]
        animals = df_s['TaxonName'].unique()
        animals = len(animals)
        count_animals.append(animals)
        # print(f'{p}[{animals}]', end=', ')
    fig = plt.figure()
    plt.bar(range(len(count_animals)), count_animals, align='center')
    firstletter_family = [f[:3]+'.' for f in data.families]
    plt.xticks(range(len(data.families)), firstletter_family, size='small')
    plt.savefig(plot_dir + 'histogram_taxfam.png', dpi=1200)
    # plt.show()


def total_average_invasion_per_year(data: InvasiveSpecies):

    avg_inv = np.zeros(data.M.shape[0])
    for idx in range(data.M.shape[0]):
        # print(np.sum(data.M[idx]))
        avg_inv[idx] = np.sum(data.M[idx])
    fig = plt.figure()
    plt.plot(data.time, avg_inv)
    plt.savefig(plot_dir + 'invasion_per_year.png', dpi=1200)
    # plt.show()

def species_region_invasion(data: InvasiveSpecies, filename='species_region_invasion.png'):
    '''Plots the amount of invasion for each species'''

    spec_inv = []
    for s in data.species:
        df_s = data.df[data.df['TaxonName'] == s]
        invasion = df_s['Region'] # get all the invasion of species s
        invasion = len(invasion)
        spec_inv.append(invasion)
    
    # Display the probability density function (pdf)
    # TODO maybe use seaborn?
    # https://stackoverflow.com/questions/71296986/how-to-draw-the-probability-density-function-pdf-plot-in-python
    fig = plt.figure()
    spec_inv.sort()
    # print(spec_inv[-3:]) # = [87, 88, 172] there's one specie that invades 172 regions
    spec_inv = spec_inv[:-1] # remove it from plot
    hmean = np.mean(spec_inv)
    hstd = np.std(spec_inv)
    pdf = stats.norm.pdf(spec_inv, hmean, hstd)
    print(f'mean: {hmean}, std{hstd}')
    plt.plot(spec_inv, pdf) 
    plt.savefig(plot_dir + filename, dpi=1200)
    # plt.show()



