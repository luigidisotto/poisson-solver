# piosson-solver

A C++ based implementation of a Poisson Equation solver performing Gauss-Seidel, Gauss-Seidel Red-Black or Jacobi methods. All the methods come with a sequential and parallel implementation.

## Requirements

* C++11 compiler, e.g GNU g++ / Intel icc. In my experiments I used the Intel compiler for perfomance reasons.

* [FastFlow](https://github.com/fastflow/fastflow), a library offering high-performance parallel patterns and building blocks in C++.

## Installation

By default, ```Makefile``` make use of GNU g++ compiler, but if you prefer using Intel icc compiler to also reproduce results, make sure to provide the name of the chosen C++11 compiler and the home directory to FastFlow, and enable the related compiler options, check CXX and OPTIMIZE_FLAGS. In that repo, find the latter as a submodule, so no worries.
Run ``` make all``` to compile the code for all the configurations configurations and version of the source code. 
Instead, if you want to build the binary for a specific method run ``` make <SOLVER>```, eg for sequential Jacobi method, do ``` make jacobi_seq```.
See below for full list of options.

### Make options

Sequential:
  * Jacobi ```jacobi_seq```
  * Gauss-Seidel ```gauss_seq```
  * Gauss-Seidel Red-Black ```gaussrb_seq```
 
Parallel:
  * Jacobi ```jacobi_par```
  * Gauss-Seidel ```gauss_par```
  * Gauss-Seidel Red-Black ```gaussrb_par```

## Basic usage

Once the code has been compiled, the generated binary code can be run provided with the following arguments ```./solver <n> <m> <eps> <nw> <sc> <bx> <by>```, where n and m represent the height and width of the grid for the discretized poisson equation, eps is the error tolerance. In the parallel versions, the option nw represents the number of threads, and sc is the FastFlow schedule type (ie data distribution policy), and bx and by are the sizes for the blocked version of the stencils, to optimize cache accesses.

To help you, the repo comes with ```run_poisson_solver.sh``` shell script, to easily build and run the desired solver.

## Theory

For further informations about basics and results, here you can find a copy of the [slides](seminar.pdf) I made for the University class on High Performance Scientific Computing @ University of Pisa.
