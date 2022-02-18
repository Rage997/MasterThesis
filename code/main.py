from combine import mesh_to_stl, populate_voxel, \
                    show_mesh, post_process
from voxel_carving import get_voxel_grid
from settings import obj_path, visualization
import logging

def main():

    # Sometimes the logger gets bugged and you need to remove all previous
    # handlers before you can log to file again. In python 3.8 they added
    # "force" param to basicConfig()
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Logging levels: CRITICAL (50), ERROR(40), WARNING(30), INFO(20), DEBUG(10), NOTSET(0)
    # If the logger value is higher than a level, no output is generated for that level
    logging.basicConfig(filename='app.log',
                    format='%(levelname)s[%(filename)s:%(lineno)d]:%(message)s',
                    filemode='w', level=logging.DEBUG, 
                    # force=True if python3.8 uncomment and remove for loop line up
                    )
    # Add stdout logging
    logging.getLogger().addHandler(logging.StreamHandler())

    logging.info('Started')   
    # meta_path = defined in settings.py
    voxel_grid, scale = get_voxel_grid(obj_path)
    voxels = voxel_grid.get_voxels()
    
    result_mesh = populate_voxel(voxels)
    result_mesh = post_process(result_mesh, scale=1)

    if visualization:
        show_mesh(result_mesh)
    mesh_to_stl(result_mesh)
    logging.info('Finished')

if __name__ == "__main__":
    main()
