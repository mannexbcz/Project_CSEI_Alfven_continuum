# Project_Alfven_continuum

## General aim of the code
This code is part of two semester projects aiming at solving the Alfvén continuum equation using a variational formulation and a Fourier series representation.
Starting from geometrically parametrized equilibria, the code has been generalized to ideal MHD axisymmetric equilibria obtained from the ideal MHD equilibrium code solver CHEASE



## Stucture of the Github
This Github project is separated into 3 Folders, two of them gathering the scripts used for the preliminary codes, and for the convergence studies & other numerical explorations written for the two projects, and a third one containing a final and functionnal version of the code.

---

### Semester 1
During this first semester, the focus was put on the resolution of the Alfvén continuum equation using a variational formulation and a Fourier series representation for simplified and geometrically parametrized equilibria. Starting from the simplest geometry with a plasma with a circular and concentric cross-section, the Shafranov shift, the elongation and the triangularity have been successively added.

The provided code include implementations of all these geometries, associated numerical computations of the Alfvén continuum as well as scripts to compare the obtained Alfvén eigenmodes and test the convergence of the different implementations.

5 versions of the code are provided, corresponding to the different geometries and an additional textbook case.
Each folder is structured in the same way and contains the following scripts
- *solver___.m* : script to use to solve the Alfvén continuum equation. Parameters can be modified directly in the file.
- *build___.m* or *matrices___.m* : functions to build the matrices of the matrix form of the Alfvén equation from the Fourier coefficients
- *eigenmodes___.m* : returns the Alfvén eigenmodes
- *get_fourier_coeff____.m* : computes the Fourier coefficients of the equilibrium coefficients

Depending on the geometry considered, some other functions may be added, allowing to compute, qbar, theta star etc...
Finally, the files named *study____.m* allow to reproduce the different convergence studies/comparison as well as all figures present in the report.

---

### Semester 2

The second semester aimed at generalizing the code to ideal MHD axisymmetric equilibria obtained from the ideal MHD equilibrium code solver CHEASE. The structure of the code and of the files is quite similar. Only one solver is however provided, the different versions of the code (Geometrical parametrization using CHEASE's coefficients, Solution using directly the equilibrium coefficients provided by CHEASE, quadratic approximation... (see report)) can be chosen while running the code. More precisely, when running the script *solver_triang.m*, the following request is printed *Method (Geometric/CHEASE/CHEASEGeom/Quadratic):* and waits for the user to specify the method. Afterwards, a window invites you if necessary to choose the *.h5* file to use.

---

### Final Code

The final code uses the solver from the second semester, with only the "CHEASE" option. Unlike the other versions, this last solver is in the form of a function that can simply be launched by calling "solver()". As before, a window then invites you to choose the *.h5* file to use.

---

__Remark :__ All functions are described in more detail directly on the files. More informations can also be found on the two reports. 
