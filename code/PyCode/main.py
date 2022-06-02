from dataset import InvasiveSpecies
import config
import make_plots
import visuals
import pickle
import os

# InvSpec = None
reset_data = True

def load_data():
    # global InvSpec
    InvSpec = InvasiveSpecies(config.filename) 
    print(f'----------Opening dataset----------------')
    InvSpec.print_info()
    # print(f'Number of species {InvSpec.Ns} and regions {InvSpec.Nr} before filtering')
    print(f'-------------Filtering elements of interest-------------')
    InvSpec.filter_data()
    InvSpec.print_info()
    # print(f'Number of species {InvSpec.Ns} and regions {InvSpec.Nr} after filtering')
    print(f'--------------------------')
    return InvSpec

def run_plots(InvSpec):
    InvSpec.build_matrix() # required for some plots
    make_plots.histogram(InvSpec)
    make_plots.total_average_invasion_per_year(InvSpec)
    make_plots.species_region_invasion(InvSpec)
    make_plots.region_species_invasion(InvSpec)


def filtering(data, n, s):
    
    # First dataset
    data.remove_irrelevant(species_tol=n, region_tol=s)
    data.print_info()

if __name__ == '__main__':

    if reset_data:
        InvSpec = load_data()
        with open('invspec_obj.pkl', 'wb+') as f:
            pickle.dump(InvSpec, f)
        # InvSpec.print_info()
    else:
        with open('invspec_obj.pkl', 'rb') as f:
            InvSpec = pickle.load(f)
    
    # Do some filtering and save depending on dataset size
    filterS = [5, 10, 15]
    filterR = [5]
    
    for s in filterS:
        for r in filterR:
            tmp = InvSpec # tbh not sure if it's deep copy
            filtering(tmp, s, r)
            # run_plots()
            tmp.build_matrix() # Need to rebuild matrix after filtering
            print(f'------ After filtering for (s,r) = ({s},{r}) ,:  ----------')
            # tmp.print_info()
            path = '../data/'+str(tmp.Ns) + '_' + str(tmp.Nr)
            if not os.path.exists(path):
                os.makedirs(path)
            name = 'matrix'
            tmp.export_matrix(path + '/' + name)
