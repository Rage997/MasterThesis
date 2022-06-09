# import cmath
from xml.dom import ValidationErr
import pandas as pd
# from sqlalchemy import column
import numpy as np
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt

path = '../data/'
# SEE: https://zenodo.org/record/4632335#.YkGzOShBz4c
filename = path+'./GlobalAlienSpeciesFirstRecordDatabase_v2.xlsx'



df_orig = pd.read_excel(filename, sheet_name=None)
# There are 3 sheets of the excel file
print(df_orig.keys())
df_orig = df_orig['GlobalAlienSpeciesFirstRecordDa']
# df_orig.info(show_counts=True)
# df_orig.head()

col = ['TaxonName', 'Family', 'LifeForm', 'Region', 'FirstRecord', 'Island']
df = df_orig[col]

df = df[(df['FirstRecord'] >= 1850) & (df['FirstRecord'] <= 2010)]
# df = df[ df['Island'].str.lower() != 'yes'] # Filter by island
# df['Island'] = 'yes'

region = df['Region'].unique()
df['Region'] = df['Region'].str.lower()
region = [reg.lower() for reg in region]

# search cross database the region
# Ref: https://www.diva-gis.org/
world_data = gpd.read_file(path+r'world.shp')

# print(world_data.head())
# world_data.plot()
# plt.show()
world_data['NAME'] = world_data['NAME'].str.lower()
world_data = world_data.replace('syrian arab republic', 'syria')
world_data = world_data.replace('libyan arab jamahiriya', 'lybia')
world_data = world_data.replace('wallis and futuna islands ', 'wallis and futuna')
world_data = world_data.replace('republic of moldova', 'moldova')
world_data = world_data.replace('united states minor outlying islands', 'us minor outlying islands')
world_data = world_data.replace('turks and caicos islands', 'turks and caicos')
world_data = world_data.replace('falkland islands (malvinas)', 'falkland islands')
world_data = world_data.replace('south georgia south sandwich islands', 'south georgia and the south sandwich islands')
world_data = world_data.replace('iran (islamic republic of)', 'iran, islamic republic of')
world_data = world_data.replace('the former yugoslav republic of macedonia', 'macedonia')
world_data = world_data.replace('palestine', 'palestine, state of')
world_data = world_data.replace("korea, democratic people's republic of", 'south korea')
world_data = world_data.replace("korea, republic of", 'north korea')
world_data = world_data.replace("viet nam", 'vietnam')
world_data = world_data.replace("democratic republic of the congo", 'congo, democratic republic of the')
world_data = world_data.replace("svalbard", 'svalbard and jan mayen')
world_data = world_data.replace("timor-leste", 'timor leste')

regions_geospatial = world_data['NAME'].unique()
# regions_geospatial = [reg.lower() for reg in regions_geospatial]

intersection = {}
intersected = []
not_intersected = []

for reg in region:
    intersection[reg] = 0
for reg2 in regions_geospatial:
    # print(reg2.lower())
    if reg2 in intersection:
        intersection[reg2] = 1
        # intersected.append(reg2)
        # print(f'Found intersection {reg2}')
    # else:
    #     print(f'Missing {reg2}')
    #     not_intersected.append(reg2)
    #     # print(f'The region {reg2} is not to be found')

# Check how many regions we found
sum = 0
for key in intersection.keys():
    sum += intersection[key]
    if intersection[key] == 1:
        intersected.append(key)
        # world_data.loc[: ,('NAME','INTERSECTED')] = 1
    else:
        # world_data.loc[: ,('NAME','INTERSECTED')] = 0
        not_intersected.append(key)

assert( ( len(intersected) + len(not_intersected) ) == len(region))
print(f'We found {sum} out of {len(region)} geospatial data')

def save_intersections():
    with open(path+"intersected.txt", "w+") as txt_file:
        for line in intersected:
            txt_file.write(line + "\n") # works with any number of elements in a line

    with open(path+"not_intersected.txt", "w+") as txt_file:
        for line in not_intersected:
            txt_file.write(line + "\n") # works with any number of elements in a line
save_intersections()

# Print geospatial data of the dataset
world_data['intersection'] = 0
world_data.loc[world_data.NAME.isin(intersected), 'intersection'] = 1

island_names = df[df['Island'] == 'yes']
island_names = island_names['Region'].unique()
island_names = [isl.lower() for isl in island_names]
# print(island_names)
with open(path+"islands.txt", "w+") as txt_file:
        for line in island_names:
            txt_file.write(line + "\n") # works with any number of elements in a line
world_data.loc[world_data.NAME.isin(island_names), 'intersection'] = 1

from matplotlib.colors import ListedColormap
road_colors = ['red', 'blue' ]

fig_path = '../../latex/figures/'
# Plot the regions that are intersected in blue and the ones mising in red
# TODO so buggy....why are there 5 bins?
fig, ax = plt.subplots(figsize  = (15, 10))
# leg_kwds_dict = {'numpoints': 2, 'labels': ['dio', 'cane']}

world_data.plot(column='intersection', 
    # cmap='tab10', 
    ax = ax,
    cmap = ListedColormap(road_colors),
    # figsize=(15, 10),
    legend=False,
    # scheme="quantiles"
    # legend_kwds=leg_kwds_dict
)
# ax.legend()
# plt.savefig(path+'dataset_intersection.png')
plt.show()

# for island
# df['Island'].replace(np.nan, 'Unknown', inplace=True)

#TODO fix database : islands and check intersections
world_data['invaded'] = 500
for region in intersected:
    invasion_count = len(df[df['Region'] == region])
    # print(f'{region} is invaded {invasion_count} times')
    world_data.loc[world_data['NAME'] == region, 'invaded']  = invasion_count 
# to_count = df_orig['Regions']
fig, ax = plt.subplots(figsize  = (15, 10))
# leg_kwds_dict = {'numpoints': 2, 'labels': ['dio', 'cane']}


most_invaded = world_data.nsmallest(20, 'invaded')
for index, row in most_invaded.iterrows():
    # print(info)
    print(row['NAME'])
    print(row['invaded'])

suisse_invasion = world_data[world_data['NAME'] == 'switzerland']
print(suisse_invasion)
# print(f'Switzerland has been invaded {suisse_invasion["invaded"]}')

world_data.plot(
    column='invaded', 
    cmap='OrRd', 
    norm=matplotlib.colors.LogNorm(vmin=world_data.invaded.min(), vmax=world_data.invaded.max()),
    ax = ax,
    legend=True,
    # scheme="quantiles"
    legend_kwds={'label': "Alien species by country from 1850 to 2010",
                        'orientation': "horizontal"}
)
# ax.legend()
plt.savefig(fig_path+'region_invasion.png')
plt.show()

#Later
#TODO plot gradient region invasion varying in time (will be used for presentation)