# N-Body Gravity Sim

This project began as an assignment for "High-Performance Computing", one of my courses at uni. The objective was to choose an algorithm that we could parallelize in one or more ways i.e. SIMD, Multiprocessing (OpenMP or MPI) or on the GPU with CUDA.

The algorithm I sought to optimize was the Barnes-Hut N-body simulation.

## How to use

### Dependencies

- g++: for compiling C++ 
- mpih: for compiling C++ with MPI using mpic++

### Building

Run `make` to compile all the binary executables.

### Running



## Slurm Scripts

I benchmarked the simulation on UQ's High-Performance Computing cluster. Jobs were run through slurm.

Note: The slurm script I have moved to the slurm_scripts folder as they are now redundant though they were intended to be run from the root directory.