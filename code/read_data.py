import pandas as pd

filename = 'GlobalAlienSpeciesFirstRecordDatabase_v2.xlsx'

pd.read_excel(filename, sheet_name=None, nrows=20)

# xl_file = pd.ExcelFile(filename)
# dfs = {sheet_name: xl_file.parse(sheet_name) 
#           for sheet_name in xl_file.sheet_names}