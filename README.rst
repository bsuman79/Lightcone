Lightcone
=========

Lightcone is a c++ library package that project the particles' positions to 
a octant spherical sky from an N-body simulation box.

Lightcone is released under the MIT software liscense (see LISCENSE).

Prerequisites
=============

Gnu Scientific Library (GSL)
A c++ compiler

Descriptions
============

The package contain 4 header files- cosmo.h, Vector.h, driver_lc.h, lightcone.h and readsnap.h. 1 input parameter file- input.param, one makefile-Makefile and one driver routine- driver_lc.cpp

1. cosmo.h- contains the class "cosmo" and the associated functions that calculates most of the cosmogical quantitites e.g., angular diamter distance and cosmic time.
3. driver_lc.h- defines and initialize some of the input parameters that 
need not be changed frequently. 
4. Makefile- makefile to generate the executable
5. input.par- file that contain the values of the input parameters.
6. driver_lc.cpp- the main driver routine that calls the class lightcone
 and project the particles from the box to the octant sphere.

 

Installation and running
========================

1. Download the package and place it anywhere you like
2. Open the Makefile, edit CFLAGS and CLIB to make sure you have the correct 
path to GSL insalled in your machine.
3. do make clean, then make lc, this creates the executable (an example is 
given in compile.sh)
4. run : ./driver_lc.exe < input.par (an example is given run_lc_omp.sh)
