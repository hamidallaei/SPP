Parallel MPI version of the self propelled particles.

To compile:
mpic++ mainsourcefile -lgsl -lcblas -O3
To run:
mpirun -np num_process a.out rho g alpha noise

The number of processes should match the input npx and npy in parameters.h file.
