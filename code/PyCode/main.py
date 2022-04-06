import dataset
import config
from REM import EKR


df = dataset.read_data(config.filename)
df = dataset.filter_data(df)


algo = EKR(df)
algo.build_matrix()
# algo.run()
