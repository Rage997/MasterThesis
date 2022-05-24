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
    print(f'Number of species {InvSpec.Ns} and regions {InvSpec.Nr} before filtering')
    InvSpec.filter_data()
    print(f'Number of species {InvSpec.Ns} and regions {InvSpec.Nr} after filtering')
    # M = InvSpec.build_matrix() # TODO only need for some plots

def run_plots():
    make_plots.histogram(InvSpec)
    make_plots.total_average_invasion_per_year(InvSpec)
    make_plots.species_region_invasion(InvSpec)
    make_plots.region_species_invasion(InvSpec)


def filtering():
    
    # First dataset
    InvSpec.remove_irrelevant(species_tol=25, region_tol=50)
    InvSpec.print_info()

if __name__ == '__main__':

    if reset_data:
        load_data()
        with open('invspec_obj.pkl', 'wb+') as f:
            pickle.dump(InvSpec, f)
    else:
        with open('invspec_obj.pkl', 'rb') as f:
            InvSpec = pickle.load(f)

        
    InvSpec.print_info()
    filtering()
    InvSpec.build_matrix() # required for some plots
    run_plots()
    InvSpec.build_matrix() # Need to rebuild matrix after filtering
    print('------ After filtering: ----------')
    InvSpec.print_info()

    InvSpec.export_matrix()
