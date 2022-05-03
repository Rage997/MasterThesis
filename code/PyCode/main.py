from dataset import InvasiveSpecies
import config
import make_plots
import visuals

InvSpec = None

def load_data():
    global InvSpec
    InvSpec = InvasiveSpecies(config.filename) 
    print(InvSpec.Ns, InvSpec.Nr)
    InvSpec.filter_data()
    InvSpec.remove_irrelevant(species_tol=3, region_tol=10)
    print(InvSpec.Ns, InvSpec.Nr)
    M = InvSpec.build_matrix() # TODO time_resolution
    
def run_plots():
    global InvSpec
    # make_plots.histogram(InvSpec)
    # make_plots.average_invasion_per_year(InvSpec)
    make_plots.species_region_invasion(InvSpec)

if __name__ == '__main__':
    load_data()
    run_plots()
