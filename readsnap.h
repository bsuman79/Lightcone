#ifndef gadget_Header_included
#define gadget_Header_included

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>


struct io_header_1
{
	int      npart[6];
	double   mass[6];
	double   time;
	double   redshift;
	int      flag_sfr;
	int      flag_feedback;
	int      npartTotal[6];
	int      flag_cooling;
	int      num_files;
	double   BoxSize;
	double   Omega0;
	double   OmegaLambda;
	double   HubbleParam; 
	char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


//int     NumPart, Ngas;


struct particle_data 
{
	float  Pos[3];
	float  Vel[3];
	float  Mass;
	int    Type;
	
	float  Rho, U, Temp, Ne;
} ;

//int *Id;

double  Time, Redshift;



/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
/*int unit_conversion(void)
{
	double GRAVITY, BOLTZMANN, PROTONMASS;
	double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
	double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
	double G, Xh, HubbleParam;
	
	int i;
	double MeanWeight, u, gamma;
	
	// physical constants in cgs units 
	GRAVITY   = 6.672e-8;
	BOLTZMANN = 1.3806e-16;
	PROTONMASS = 1.6726e-24;
	
	// internal unit system of the code 
	UnitLength_in_cm= 3.085678e21;   //  code length unit in cm/h 
	UnitMass_in_g= 1.989e43;         //  code mass unit in g/h 
	UnitVelocity_in_cm_per_s= 1.0e5;
	
	UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
	UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
	UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
	UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
	
	G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);
	
	
	Xh= 0.76;  // mass fraction of hydrogen 
	HubbleParam= 0.65;
	
	
	for(i=1; i<=NumPart; i++)
    {
		if(P[i].Type==0)  // gas particle 
		{
			MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;
			
			// convert internal energy to cgs units 
			
			u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;
			
			gamma= 5.0/3;
			
			// get temperature in Kelvin 
			
			P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
		}
    }
}
*/
/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(int NumPart, particle_data *P)
{
	printf("allocating memory...\n");
	
	//if(!(P=malloc(NumPart*sizeof(struct particle_data))))
      if(!(P=new struct particle_data[NumPart]))
    {
		fprintf(stderr,"failed to allocate memory.\n");
		exit(0);
    }
	
	P--;   /* start with offset 1 */
	
	
	//if(!(Id=malloc(NumPart*sizeof(int))))
      /*if(!(Id= new int [NumPart]))
    {
		fprintf(stderr,"failed to allocate memory.\n");
		exit(0);
    }
	*/
	//Id--;   /* start with offset 1 */
	
	//printf("allocating memory...done\n");
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int get_partnum(char *fname)
{
       FILE *fd;
       char   buf[200];
       int    i,j,k,dummy,ntot_withmasses;
       int NumPart, Ngas;   
       sprintf(buf,"%s",fname);
       if(!(fd=fopen(buf,"r")))
	{
		printf("can't open file `%s`\n",buf);
		exit(0);
	}
	
	printf("reading `%s' ...\n",buf); fflush(stdout);
	fread(&dummy, sizeof(dummy), 1, fd);
	fread(&header1, sizeof(header1), 1, fd);
	fread(&dummy, sizeof(dummy), 1, fd);
	
	for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	Ngas= header1.npart[0];
	
	for(k=0, ntot_withmasses=0; k<5; k++)
	{
		if(header1.mass[k]==0)
			ntot_withmasses+= header1.npart[k];
	}
    fclose(fd);
    return NumPart;
	
}

int load_snapshot(char *fname,particle_data *P)
{
	FILE *fd;
	char   buf[200];
	int    i,j,k,dummy,ntot_withmasses;
	int    t,n,off,pc,pc_new,pc_sph;
	int NumPart, Ngas;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
	float min[3], max[3];
	min[0]=1e6; min[1]= 1e6; min[2]= 1e6;
	max[0]=-1.0; max[1]= -1.0; max[2]= -1.0;
	sprintf(buf,"%s",fname);
	
	if(!(fd=fopen(buf,"r")))
	{
		printf("can't open file `%s`\n",buf);
		exit(0);
	}
	
	printf("reading `%s' ...\n",buf); fflush(stdout);
	
	fread(&dummy, sizeof(dummy), 1, fd);
	fread(&header1, sizeof(header1), 1, fd);
	fread(&dummy, sizeof(dummy), 1, fd);
	
	for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	Ngas= header1.npart[0];
	
	for(k=0, ntot_withmasses=0; k<5; k++)
	{
		if(header1.mass[k]==0)
			ntot_withmasses+= header1.npart[k];
	}
	 
	//allocate_memory();
	
	SKIP;
	for(k=0, pc_new=1;k<6;k++)
	{       // printf("%d\n",header1.npart[k]);
		for(n=0;n<header1.npart[k];n++)
	    {
			fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
			//if(n<50) printf("%f %f %f\n",P[pc_new].Pos[0], P[pc_new].Pos[1], P[pc_new].Pos[2]);	
			if(min[0] > P[pc_new].Pos[0]) min[0]=P[pc_new].Pos[0];
			if(max[0] < P[pc_new].Pos[0]) max[0]=P[pc_new].Pos[0];
			if(min[1] > P[pc_new].Pos[0]) min[1]=P[pc_new].Pos[1];
			if(max[1] < P[pc_new].Pos[0]) max[1]=P[pc_new].Pos[1];
			if(min[2] > P[pc_new].Pos[0]) min[2]=P[pc_new].Pos[2];
			if(max[2] < P[pc_new].Pos[0]) max[2]=P[pc_new].Pos[2];
			pc_new++;
	    }
	}
	SKIP;
	
	SKIP;
	for(k=0, pc_new=1;k<6;k++)
	{
	        for(n=0;n<header1.npart[k];n++)
			fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
                pc_new++;
	}
	
	fclose(fd);
    
	//printf("%d\n",NumPart);
	printf("%f %f %f %f %f %f\n", min[0], max[0], min[1], max[1], min[2], max[2]);
	Time= header1.time;
	Redshift= header1.time;
}


void deletepos(particle_data *P)
{
 delete[] P;
 //delete[] Id;
}





/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 *//*
int reordering(void)
{
	int i,j;
	int idsource, idsave, dest;
	struct particle_data psave, psource;
	
	
	printf("reordering....\n");
	
	for(i=1; i<=NumPart; i++)
    {
		if(Id[i] != i)
		{
			psource= P[i];
			idsource=Id[i];
			dest=Id[i];
			
			do
			{
				psave= P[dest];
				idsave=Id[dest];
				
				P[dest]= psource;
				Id[dest]= idsource;
				
				if(dest == i) 
					break;
				
				psource= psave;
				idsource=idsave;
				
				dest=idsource;
			}
			while(1);
		}
    }
	
	printf("done.\n");
	
	Id++;   
	free(Id);
	
	printf("space for particle ID freed\n");
}*/
#endif
