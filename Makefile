# mac
#CC = g++ # uncomment this line for carver, comment out for hopper
CFLAGS = $(GSL) -fopenmp -O3  -ffast-math -funroll-loops   #-I/usr/local/include -fopenmp
#CFLAGS = $(GSL) -O3  -ffast-math -funroll-loops 
CLIBS = $(GSL) #-L/usr/local/lib -lgsl -lgslcblas  -lm

.cpp.o :
	$(CC) $(CFLAGS)  -c $<


lc: driver_lc.o
	$(CC) -o driver_lc.exe driver_lc.o $(CLIBS)

clean : 
	/bin/rm -f *.o *.exe core.*

