import os

# Voxelization settings
data_folder = '../data/'
output_filename = os.path.abspath(f'{data_folder}result.stl')
camera_path = os.path.abspath(f'{data_folder}sphere.ply')
meta_path =  os.path.abspath(f'{data_folder}meta.stl')
obj_path = os.path.abspath(f'{data_folder}cube.stl')

visualization = False
cubic_size = 1.0 # controls voxel grid size
voxel_resolution = 1
voxel_size = cubic_size / voxel_resolution 

voxel_per_dimension = 2*cubic_size/voxel_size
total_voxel = voxel_per_dimension**3
# voxel per dimension = 1/
# Final object size
# size = 
# TODO voxel size = metastrucutre size
