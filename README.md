# Master Thesis


Repository structure:

1. Metastructure: this folder contains everything related to the generation of the metastructure.
2. Macrostructure: this folder contains everything related to the assembly of the final structure (i.e. by composing the metastructure)
3. Testing: some jupyter notebooks to test functions and libraries
4. Data: all the generated 3D files

Todos:

[ ] Optimise the number of vertices metastrucuture

[ ] Logging: number of voxels, time to voxelise, to generate metastructure etc

[ ] Multithreading final mesh generation

[ ] Metastructure on boundaries?

# Installation
I used conda enviroment to manage packages. There are two environments, one for the generate the metastructure and one to generate the final mesh.

You can install them both by running:

```conda create --name cadquery --file ./metastructure/env.txt```
```conda create --name open3d --file ./code/env.txt```

# Performance

The code is fast. To fill a voxelGrid of 715597 voxels, it took 43 seconds on my machine (Ryzen 3700).

# Printing

Printing the metastructure is exceptionally hard due to its fine details. Lower your wall line width to something like half the wall thickness (i.e. 0.225 mm) and see if that works. With a standard 0.4 mm nozzle I've had success printing tiny details with 0.2 mm line width or smaller. In Cura (not sure about other slicers) you also need to set the "Outer wall inset" to zero.
You also need to drastically decrease the print speed. I've been printing at ~12 mm/s. After the first layer, you can slightly increase the speed.
You may also 