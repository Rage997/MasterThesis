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
    print(InvSpec.Ns, InvSpec.Nr)
    InvSpec.filter_data()
    # InvSpec.remove_irrelevant(species_tol=3, region_tol=10)
    print(InvSpec.Ns, InvSpec.Nr)
    # M = InvSpec.build_matrix() # TODO only need for some plots

def run_plots():
    make_plots.histogram(InvSpec)
    make_plots.total_average_invasion_per_year(InvSpec)
    make_plots.species_region_invasion(InvSpec)

def filtering():
    InvSpec.remove_irrelevant(species_tol=25, region_tol=10)
    InvSpec.print_info()
    make_plots.species_region_invasion(InvSpec, filename='species_region_invasion_filtering.png')

if __name__ == '__main__':

    if reset_data:
        load_data()
        with open('invspec_obj.pkl', 'wb+') as f:
            pickle.dump(InvSpec, f)
    else:
        with open('invspec_obj.pkl', 'rb') as f:
            InvSpec = pickle.load(f)

        
    InvSpec.print_info()
    # run_plots()
    filtering()
    InvSpec.build_matrix() # Need to rebuild matrix after filtering
    print('------ After filtering: ----------')
    InvSpec.print_info()

    # InvSpec.export_matrix()
