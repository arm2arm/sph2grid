Install:
./configure
To use enable SMP parallel version of the code set :
setenv CXXFLAGS " -openmp -D_OPENMP "
setenv CXX icpc
Then run ./configure
Usage:
SPH2GRID project files
    In: gadget file [Fname]
        grid size for 3 [NX NY NZ]
        center of the grid  [XC YC ZC]
        Radius of the center [RC]
    OUT: Gridfile [Fname].grid
      
The format of the *.grid file is following:
Fortran 77 format similar like Gadget-2
HEAD : size 256 ( 3 INT[NX NY NZ], 3 FLOAT [XC YC ZC], 1 FLOAT[Rc])
GRID : FLOAT Arr3D [Nx][Ny][Nz]
/////////////////////////////////////////////////////////////////////////////
