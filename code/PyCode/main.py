from dataset import InvasiveSpecies
import config
from REM import EKR


InvSpec = InvasiveSpecies(config.filename) 
# InvSpec.filter_data()
M = InvSpec.build_matrix()
InvSpec.export_to_R()

algo = EKR(M)
# algo.build_matrix()
# algo.run()
