salloc -J masterthesis -t 600 --partition=gpu
module load gcc/10.1.0
module load cudatoolkit
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/
