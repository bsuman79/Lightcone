#!/bin/tcsh
module swap PrgEnv-pgi PrgEnv-gnu
module load gsl
make clean
make lc
