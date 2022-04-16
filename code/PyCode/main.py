from dataset import InvasiveSpecies
import config
from REM import EKR

InvSpec = InvasiveSpecies(config.filename) 
print(InvSpec.Ns, InvSpec.Nr)
InvSpec.filter_data()
InvSpec.remove_irrelevant(species_tol=3, region_tol=10)
print(InvSpec.Ns, InvSpec.Nr)
M = InvSpec.build_matrix(time_resolution=10)
# # InvSpec.export_to_R()

algo = EKR(M)
algo.run(100)
