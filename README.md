# Project_CSEI_Alfven_continuum

## General aim of the code
This code is part of a semester project aiming at solving the Alfvén continuum equation using a variational formulation and a Fourier series representation. Starting from the simplest geometry with a plasma with a circular and concentric cross-section, the Shafranov shift, the elongation and the triangularity have been successively added.

The provided code include implementations of all these geometries, associated numerical computations of the Alfvén continuum as well as scripts to compare the obtained Alfvén eigenmodes and test the convergence of the different implementations.

## Structure of the code
5 versions of the code are provided, corresponding to the different geometries and an additional textbook case.
Each folder is structured in the same way and contains the following scripts
- *solver___.m* : script to use to solve the Alfvén continuum equation. Parameters can be modified directly in the file.
- *build___.m* or *matrices___.m* : functions to build the matrices of the matrix form of the Alfvén equation from the Fourier coefficients
- *eigenmodes___.m* : returns the Alfvén eigenmodes
- *get_fourier_coeff____.m* : computes the Fourier coefficients of the equilibrium coefficients

Depending on the geometry considered, some other functions may be added, allowing to compute, qbar, theta star etc...
Finally, the files named *study____.m* allow to reproduce the different convergence studies/comparison as well as all figures present in the report.

All functions are described in more detail directly on the files. More informations can be found on the project report. Note also that the equations referenced in the documentation of the code refer to the equations derived in the report.

