SPP
===

Serial code plus an MPI main that runs multiple serial programs with different parameters.

To compile the simulator (single):
g++ -lgsl -lcblas -O3 main.cpp

To compile the simulator (Multiple):
mpic++ -lgsl -lcblas -O3 mpi-main.cpp
