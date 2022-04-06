from dataset import InvasiveSpecies
import config
from REM import EKR


InvSpec = InvasiveSpecies(config.filename) 
# InvSpec.filter_data()
M = InvSpec.build_matrix()


# algo = EKR(df)
# algo.build_matrix()
# algo.run()
