#include <stdio.h>
#include <math.h>
#include "lightcone.h"

//#include "BasicDefinition.h"
#include "readsnap.h"
#include "driver_lc.h"

#define rhoc 2.77536627e11
#define  pi 3.141593
#define clight 3.0e5
using namespace std;
int main(int argc, char** argv)
{
  /*****************************/
  /*cosmology params*/
 /*****************************/

        cout<< "read in parameters"<<endl;
        cin >> Ob;
        cin >> Om;
        cin >> wt;
        cin >> H0;
        cin >> boxsize;
        cin >> numsnaps;
        cin >> numfile;
        cin >> np;
        cin >> snapindx;
        cin >> outFile;
        cin >> inFile;
        sprintf(path, inFile);//"/project/projectdirs/hacc/Coyote/M000/G004/Gadget_");


  //sprintf(path, "/project/projectdirs/hacc/Coyote/M000/G003/Gadget_");

  //string outFile = "/global/project/projectdirs/hacc/lightcone/jnk_";
  
  lightcone lc(H0, Om, Ob, 0.0, wt, snapindx, scalefac[snapindx-1], scalefac[snapindx], scalefac[snapindx+1], boxsize);
   cout<<"this snapshot: "<<1/(scalefac[snapindx])-1<<" "<< lc.get_comovinglow()<<" "<<lc.get_comovinghigh()<<" "<<lc.get_comovingav()<<" "<<lc.get_comovingwidth()<<endl;

   double mp= pow(boxsize/np,3)*Om*rhoc*H0*H0/1.0e4;
   cout<<"mass of a particle="<<mp<<endl;
   int j, num=0;
   int numpart;//= floor(2*1e9/numfile); 
    particle_data *P;
  for(j=0;j<numfile;j++){
    sname<<outFile<<addtofile[snapindx]<<"."<<j<<".lc";
    ofstream outfile(sname.str().c_str(), ios::out|ios::binary);  
    sprintf(input_fname, "%s%2.4f.%d", path, scalefac[snapindx], j);
    numpart= get_partnum(input_fname);
    cout<<numpart<<endl;
    if (j==0){
       if(!(P=new struct particle_data[numpart+numpart]))
          {
                fprintf(stderr,"failed to allocate memory.\n");
                exit(0);
          }
     
     P--; 
     }
    // allocate_memory(numpart, P);
    load_snapshot(input_fname, P);

    particleCount= numpart;
    cout<<"particle # in file "<<j<<"="<<particleCount<<endl;
    particle_lc_file=0;
    for(long ii=0;ii<particleCount;ii++) {
        pos[0]= (double)P[ii+1].Pos[0]*100.0/H0;
        pos[1]= (double)P[ii+1].Pos[1]*100.0/H0;
        pos[2]= (double)P[ii+1].Pos[2]*100.0/H0;
        //cout<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
        lc.map_a_particle_position_to_snapshot_slice(pos,outfile, num);        
        particle_lc_file += num;
     }
    particle_lc += particle_lc_file;
    cout<<"particles in lightcone in file "<<j<<"="<<particle_lc_file<<endl;
    lc.close_output_file(outfile);
    sname.str("");
    //delete[] P;
    //deletepos(P);
   cout<<"done file "<<j<<endl;
  }
   
  double total_mass_lc=particle_lc*mp;
  double volume_lc= (pi/6.0*(pow(lc.get_comovinghigh(),3)-pow(lc.get_comovinglow(),3)));
  double density_lc= total_mass_lc/volume_lc;
  cout<<"total particles in lightcone="<<particle_lc<<" density= "<<density_lc<< " Msun/Mpc^3"<<endl;
  cout<<" ratio wrt to the mean density="<< density_lc/(Om*rhoc*H0*H0/1.0e4)<<endl;
  cout<<" volume of the shell="<<volume_lc*1e-9<<" Gpc^3"<<endl;
  return 1;
}
