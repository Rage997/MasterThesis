from dataset import InvasiveSpecies
import config
import make_plots
import visuals
import pickle

InvSpec = None
reset_data = True

def load_data():
    global InvSpec
    InvSpec = InvasiveSpecies(config.filename) 
    print(f'----------Opening dataset----------------')
    InvSpec.print_info()
    # print(f'Number of species {InvSpec.Ns} and regions {InvSpec.Nr} before filtering')
    print(f'-------------Filtering elements of interest-------------')
    InvSpec.filter_data()
    InvSpec.print_info()
    # print(f'Number of species {InvSpec.Ns} and regions {InvSpec.Nr} after filtering')
    print(f'--------------------------')

def run_plots():
    InvSpec.build_matrix() # required for some plots
    make_plots.histogram(InvSpec)
    make_plots.total_average_invasion_per_year(InvSpec)
    make_plots.species_region_invasion(InvSpec)
    make_plots.region_species_invasion(InvSpec)


def filtering(n, s):
    
    # First dataset
    InvSpec.remove_irrelevant(species_tol=n, region_tol=s)
    InvSpec.print_info()

if __name__ == '__main__':

    if reset_data:
        load_data()
        with open('invspec_obj.pkl', 'wb+') as f:
            pickle.dump(InvSpec, f)
    else:
        with open('invspec_obj.pkl', 'rb') as f:
            InvSpec = pickle.load(f)

        
    filtering(15, 15)
    # run_plots()
    InvSpec.build_matrix() # Need to rebuild matrix after filtering
    print('------ After filtering: ----------')
    InvSpec.print_info()

    InvSpec.export_matrix()
