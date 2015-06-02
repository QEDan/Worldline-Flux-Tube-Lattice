# Worldline-Flux-Tube-Lattice
C and CUDA code for computing the effective actions of cylindrically symmetric magnetic field configurations using worldline numerics parallelized on Nvidia GPUs.

## Documentation:
The most detailed documentation for this code is provided in my PhD thesis:
http://arxiv.org/abs/1209.4409

Note that this code is functional, but not mature. It was developed primarily in 2011 as part of a physics PhD project. The following important improvements are still outstanding:
- Use of CUDA memory heirarchy should be improved. The worldline data did not fit into shared memory on the target devices, but on newer devices this should not be a problem. Occupancy can likely be improved by moving numerical constants from registers into constant memory.
- Asynchronous CPU and GPU computing and memory transfers can be improved.
- CPU-side integrations can be improved by using optimized third-party libraries. There may also be optimized integration libraries suitable for use in the GPU kernels.

## Installation:
The code uses few dependencies. There is a basic Makefile that has worked on at least two different Linux systems. 

## Citations:
Casimir interactions between magnetic flux tubes in a dense lattice

Dan Mazur and Jeremy S. Heyl

Phys. Rev. D 91, 065019 – Published 16 March 2015

http://journals.aps.org/prd/abstract/10.1103/PhysRevD.91.065019

