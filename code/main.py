from combine import mesh_to_stl, populate_voxel, \
                    show_mesh, post_process
from voxel_carving import get_voxel_grid
from settings import obj_path


def main():
    # meta_path = defined in settings.py
    voxel_grid, scale = get_voxel_grid(obj_path)
    voxels = voxel_grid.get_voxels()
    
    result_mesh = populate_voxel(voxels)
    result_mesh = post_process(result_mesh, scale=1)

    show_mesh(result_mesh)
    mesh_to_stl(result_mesh)

if __name__ == "__main__":
    main()
