#PBS -A hacc
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=11:00:00
#PBS -N lc
#PBS -e lc.z4.err
#PBS -o lc.z4.out
#PBS -V

cd $PBS_O_WORKDIR

setenv OMP_NUM_THREADS 1
#aprun -n 1 -N 1 -d 16 
./driver_lc_openmp.exe 10
