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