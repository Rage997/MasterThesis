import os

# Voxelization settings
data_folder = '../data/'
output_filename = os.path.abspath(f'{data_folder}result.stl')
camera_path = os.path.abspath(f'{data_folder}sphere.ply')
meta_path =  os.path.abspath(f'{data_folder}meta.stl')
obj_path = os.path.abspath(f'{data_folder}suzanne.stl')

visualization = True
cubic_size = 2.0
voxel_resolution = 2**6
