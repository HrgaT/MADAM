unset LD_PRELOAD
export LD_LIBRARY_PATH=${PWD}/SuiteSparse/lib:${LD_LIBRARY_PATH}
export SUITESPARSE=${PWD}/SuiteSparse
module load OpenBLAS 
module load OpenMPI/3.1.4-GCC-8.3.0
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
