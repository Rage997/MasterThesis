import numpy as np
from voxel_carving import get_voxel_grid
import open3d as o3d
# import threading
import timeit
from settings import output_filename, meta_path, obj_path


def populate_voxel(voxel):
    start = timeit.default_timer()
    result_mesh = meta_to_voxel(voxel)
    stop = timeit.default_timer()
    print('Generating the mesh took: ', stop - start)
    return result_mesh

def meta_to_voxel(voxel: o3d.geometry.Voxel):
    '''Insert a metastructure in a voxel'''
    
    result_mesh = o3d.geometry.TriangleMesh()
    for v in voxel:
        meta_obj = o3d.io.read_triangle_mesh(meta_path)
        # TODO do not preprocess?
        # meta_obj =  preprocess(meta_obj)
        # TODO we can use the mesh color later for FEM (?)
        # meta_obj.paint_uniform_color(v.color)
        meta_obj.translate(v.grid_index, relative=False)
        result_mesh += meta_obj

    return result_mesh

def post_process(obj: o3d.geometry.TriangleMesh, scale = 1):
    obj.translate([0, 0, 0], relative=True)
    obj.scale(scale, center=[0, 0, 0])
    # obj.translate(voxel_grid.origin, relative=True)
    obj.merge_close_vertices(0.0000001)
    obj = o3d.geometry.TriangleMesh.compute_triangle_normals(obj)
    return obj

def mesh_to_stl(mesh, path = output_filename):
    o3d.io.write_triangle_mesh(path, mesh)

def show_mesh(mesh):
    o3d.visualization.draw_geometries([mesh])

if __name__ == '__main__':
    # For testing

    voxel_grid = get_voxel_grid(obj_path)
    voxels = voxel_grid.get_voxels()

    result_mesh = o3d.geometry.TriangleMesh()

    start = timeit.default_timer()
    meta_to_voxel(voxels)
    stop = timeit.default_timer()
    print('Generating the mesh took: ', stop - start)

    result_mesh = post_process(result_mesh)
    o3d.visualization.draw_geometries([result_mesh])
    # o3d.io.write_triangle_mesh("../data/"+"combined.stl", result_mesh)