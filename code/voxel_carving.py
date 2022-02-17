import open3d as o3d
import numpy as np
# import copy
# import os
from settings import output_filename, camera_path, visualization, cubic_size, voxel_resolution

def xyz_spherical(xyz):
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    r = np.sqrt(x * x + y * y + z * z)
    r_x = np.arccos(y / r)
    r_y = np.arctan2(z, x)
    return [r, r_x, r_y]

def get_rotation_matrix(r_x, r_y):
    rot_x = np.asarray([[1, 0, 0], [0, np.cos(r_x), -np.sin(r_x)],
                        [0, np.sin(r_x), np.cos(r_x)]])
    rot_y = np.asarray([[np.cos(r_y), 0, np.sin(r_y)], [0, 1, 0],
                        [-np.sin(r_y), 0, np.cos(r_y)]])
    return rot_y.dot(rot_x)

def get_extrinsic(xyz):
    rvec = xyz_spherical(xyz)
    r = get_rotation_matrix(rvec[1], rvec[2])
    t = np.asarray([0, 0, 2]).transpose()
    trans = np.eye(4)
    trans[:3, :3] = r
    trans[:3, 3] = t
    return trans

def get_scale(model: o3d.geometry.TriangleMesh):
    min_bound = model.get_min_bound()
    max_bound = model.get_max_bound()
    scale = np.linalg.norm(max_bound - min_bound) / 2.0
    return scale

def preprocess(model: o3d.geometry.TriangleMesh) -> o3d.geometry.TriangleMesh:
    """Pre-process an object by centering it to [0, 0, 0] and uniformly scaling it to [1, 1, 1]
     
    Args:
        model (o3d.geometry.TriangleMesh): The model to pre-process

    Returns:
        o3d.geometry.TriangleMesh: Pre-processed model
    """
    global scale

    min_bound = model.get_min_bound()
    max_bound = model.get_max_bound()
    center = min_bound + (max_bound - min_bound) / 2.0
    scale = np.linalg.norm(max_bound - min_bound) / 2.0
    vertices = np.asarray(model.vertices)
    vertices -= center
    model.vertices = o3d.utility.Vector3dVector(vertices / scale)

    # save for debugging
    # o3d.io.write_triangle_mesh(f'../data/preprocess.stl', model)
    return model

def voxel_carving(obj: o3d.geometry.TriangleMesh,
                  camera_path,
                  cubic_size,
                  voxel_resolution
                  ):
    
    obj.compute_vertex_normals()
    camera_sphere = o3d.io.read_triangle_mesh(camera_path)

    # Setup dense voxel grid
    voxel_carving = o3d.geometry.VoxelGrid.create_dense(
        width=cubic_size,
        height=cubic_size,
        depth=cubic_size,
        voxel_size=cubic_size / voxel_resolution,
        origin=[-cubic_size / 2.0, -cubic_size / 2.0, -cubic_size / 2.0],
        color=[1.0, 0.7, 0.0])

    # Pre-process geometry
    camera_sphere = preprocess(camera_sphere)
    obj = preprocess(obj)

    # setup visualizer to render depthmaps
    vis = o3d.visualization.Visualizer()
    vis.create_window(width=300, height=300, visible=False)
    vis.add_geometry(obj)
    vis.get_render_option().mesh_show_back_face = True
    ctr = vis.get_view_control()
    param = ctr.convert_to_pinhole_camera_parameters()

    # carve voxel grid
    pcd_agg = o3d.geometry.PointCloud()
    centers_pts = np.zeros((len(camera_sphere.vertices), 3))
    for cid, xyz in enumerate(camera_sphere.vertices):
        # get new camera pose
        trans = get_extrinsic(xyz)
        param.extrinsic = trans
        c = np.linalg.inv(trans).dot(np.asarray([0, 0, 0, 1]).transpose())
        centers_pts[cid, :] = c[:3]
        ctr.convert_from_pinhole_camera_parameters(param)

        # capture depth image and make a point cloud
        vis.poll_events()
        vis.update_renderer()
        depth = vis.capture_depth_float_buffer(False)
        pcd_agg += o3d.geometry.PointCloud.create_from_depth_image(
            o3d.geometry.Image(depth),
            param.intrinsic,
            param.extrinsic,
            depth_scale=1)

        voxel_carving.carve_depth_map(o3d.geometry.Image(depth), param)
        if cid%100 == 0:
            print("Carve view %03d/%03d" % (cid + 1, len(camera_sphere.vertices)))
    vis.destroy_window()

    voxel_surface = o3d.geometry.VoxelGrid.create_from_point_cloud_within_bounds(
        pcd_agg,
        voxel_size=cubic_size / voxel_resolution,
        min_bound=(-cubic_size / 2, -cubic_size / 2, -cubic_size / 2),
        max_bound=(cubic_size / 2, cubic_size / 2, cubic_size / 2))
    voxel_carving_surface = voxel_surface + voxel_carving

    return voxel_carving_surface, voxel_carving, voxel_surface

def get_voxel_grid(obj_path: str):
    """Given a path to a stl, it returns a voxel grid in unit coordinates and the
        world matrix of the original mesh

    Args:
        obj_path (_type_): path to stl

    Returns:
        o3d.VoxelGrid: VoxelGrid in unit coordinates

    """    
    
    global cubic_size, voxel_resolution, camera_path, scale

    base_obj = o3d.io.read_triangle_mesh(obj_path)
    voxel_grid, _, _ = voxel_carving(
        base_obj, camera_path, cubic_size, voxel_resolution)

    print(f'The scale of the object is: {scale}')
    return voxel_grid, scale

if __name__ == '__main__':
    # Testing library

    voxel_grid = get_voxel_grid('../data/suzanne.stl')

    print("combined voxels (carved + surface)")
    print(voxel_grid)
    o3d.visualization.draw_geometries([voxel_grid])