from combine import mesh_to_stl, populate_voxel, show_mesh
from voxel_carving import get_voxel_grid
from settings import obj_path


def main():
    # meta_path = defined in settings.py
    voxel_grid = get_voxel_grid(obj_path)
    voxels = voxel_grid.get_voxels()
    
    result_mesh = populate_voxel(voxels)

    show_mesh(result_mesh)
    mesh_to_stl(result_mesh)

if __name__ == "__main__":
    main()
