#ifndef driverlc_Header_included
#define driverlc_Header_included


/*cosmology params*/
float Ob, Om, wt, H0; // model 0
// float Ob=0.0647;float Om=0.431;float wt=-0.82; float H0=59.77;//model 1
float boxsize;
int numsnaps;
int numfile;
int np; //cube of this num= particles in sim

float scalefac[]={0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3333, 0.2857, 0.25, 0.2, 0.13}; //snapshots scale fact  of each gadget run
string addtofile[]={"1.0000","0.9000", "0.8000", "0.7000", "0.6000","0.5000","0.4000","0.3333","0.2857","0.2500","0.2000"};
/*****************************/
int snapindx; //which z, need to be >0 
ostringstream sname;
string outFile;
char inFile[200];
char path[200]/*, basename[200]*/, input_fname[200];
float pos[3];
long particle_lc=0, particle_lc_file;
long particleCount;
#endif
