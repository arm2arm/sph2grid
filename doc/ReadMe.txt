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
    or
    OUT: sph2grid.png  
The format of the *.grid file is following:
Fortran 77 format similar like Gadget-2
HEAD : size 256 ( 3 INT[NX NY NZ], 3 FLOAT [XC YC ZC], 1 FLOAT[Rc])
GRID : FLOAT Arr3D [Nx][Ny][Nz]
/////////////////////////////////////////////////////////////////////////////
The code is recognising the following environement variables:
SPH2GRID_DISPLAY 0/1 - hide/show the interactive image window
SPH2GRID_LUT     1-43 - select the color table, default=14
SPH2GRID_XYZ    0/1  - output the XY and XZ projections default=1
SPH2GRID_MINRHO      - set the minimum density plot in the logarithmic units default=-7.0f
SPH2GRID_MAXRHO      - set the maximum density plot in the logarithmic units default=3.0f

rendering style:
SPH2GRID_MAXINTENSITY 0/1 - plot the  particles sorted by density, default=0




