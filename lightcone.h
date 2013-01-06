#ifndef lightcone_Header_included
#define lightcone_Header_included

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream> 
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "cosmo.h"
#include "Vector.h"

using namespace std;

class lightcone{
 protected:
  float H0, Omega_M, Omega_b, rho_crit, a, Omega_k, Omega_L, wt;
  float PI, m_sun, G, mpc, gnewt;
  float scalefac_low, scalefac, scalefac_high, delta_scalefac_low, delta_scalefac_high;
  double comov_low, comov_high, comov_av, comov_width;
  //float pos[3];
  float boxsize;
 
  //vector<float> posx, posy, posz;
 
  void comovingdist(cosmo & cosm_model) {
       delta_scalefac_low= 0.5*(scalefac_low-scalefac);
       delta_scalefac_high= 0.5*(scalefac-scalefac_high);
       comov_low=cosm_model.ang_diam(1.0/(scalefac+delta_scalefac_low)-1.0)/(scalefac+delta_scalefac_low);
       comov_high=cosm_model.ang_diam(1.0/(scalefac-delta_scalefac_high)-1.0)/(scalefac-delta_scalefac_high);
       comov_av=0.5*(comov_low+comov_high);
       comov_width= comov_high-comov_low;
    }
  
  int num_boxes(){
     return (int) ceil(comov_high/boxsize);
   }

 public:
 lightcone(float inp1, float inp2, float inp3, float inp4, float inp5, int snapindx, float inp6, float inp7, float inp8, float inp9)
 {
    // cosmological params at z =0
     H0 = inp1;
     Omega_M = inp2;
     Omega_b = inp3;
     Omega_k = inp4;
     Omega_L = 1.0- Omega_M- Omega_k;
     wt= inp5;
     scalefac_low= inp6;
     scalefac= inp7;
     scalefac_high= inp8;
     boxsize= inp9;
     cosmo cosm_model(H0, Omega_M, 0, 0, wt);
     comovingdist(cosm_model);
     PI = 4.*atan(1.);
     m_sun = 1.98892e30;
     G = 6.67e-11/1.0e9; // in km^3 kg^-1 s^-2
     mpc = 3.0857e19; // in km
     gnewt = G*m_sun/mpc; //for cosm params (in km^2 Mpc msun^-1 s^-2)
  }

 /* void read_part_pos(float *inp) {
     pos[0]= inp[0];
     pos[1]= inp[1];
     pos[2]= inp[2];
     }*/
  
  void map_a_particle_position_to_stacked_boxes(float pos[3], vector<float> &xx, vector<float> &yy, vector<float> &zz){
    CVector C;
    int ctr=0;
    int N= num_boxes();
    float tmp1, tmp2, tmp3;
    float radius;
    // cout<<"number of boxes= "<<N<<endl;
    for(int i=0; i< N; i++){
      for(int j=0; j< N; j++){
        for(int k=0; k< N; k++){
	  CVector A(i,j,k, boxsize);
          CVector B(pos[0], pos[1], pos[2],1);
          C= A + B;
          tmp1= C.x;
          tmp2= C.y;
          tmp3= C.z;
          radius= sqrt(pow(tmp1,2) + pow(tmp2,2)+ pow(tmp3,2));
	 
          if(radius>= comov_low && radius <comov_high) {
	     xx.push_back(tmp1); yy.push_back(tmp2); zz.push_back(tmp3); 
	     //cout<<" radius= "<<radius<<" "<<tmp1<<" "<<tmp2<<" "<<tmp3<<endl;
           ctr++;
          }  
          //cout <<" "<< posx[ctr]<<" "<<posy[ctr]<<" "<<posz[ctr]<<endl;
         
         }
       }
    }
    cout<<"num= "<<ctr<<endl;
    xx.resize(ctr); yy.resize(ctr); zz.resize(ctr);  
  }
 

  void map_a_particle_position_to_snapshot_slice(float pos[3], ofstream &outfile, int &num){
    CVector C;
    int ctr=0;
    int N= num_boxes();
    float tmp1, tmp2, tmp3;
    float radius;
    // cout<<"number of boxes= "<<N<<endl;
    for(int i=0; i< N; i++){
      for(int j=0; j< N; j++){
        for(int k=0; k< N; k++){
	  CVector A(i,j,k, boxsize);
          CVector B(pos[0], pos[1], pos[2],1);
          C= A + B;
          tmp1= C.x;
          tmp2= C.y;
          tmp3= C.z;
          radius= sqrt(pow(tmp1,2) + pow(tmp2,2)+ pow(tmp3,2));
          if(radius>= comov_low && radius <comov_high) {
              outfile.write(reinterpret_cast<const char*>(&tmp1), sizeof (float));
              outfile.write(reinterpret_cast<const char*>(&tmp2), sizeof (float));
              outfile.write(reinterpret_cast<const char*>(&tmp3), sizeof (float)); 
	      // cout<<" radius= "<<radius<<" "<<tmp1<<" "<<tmp2<<" "<<tmp3<<endl;
           ctr++;
          }  
          //cout <<" "<< posx[ctr]<<" "<<posy[ctr]<<" "<<posz[ctr]<<endl;
         
         }
       }
    }
    //cout<<"num= "<<ctr<<endl;
    num=ctr;
  }
 

  void mapped_particles_in_snapshot_slice( vector<float> &posx, vector<float> &posy, vector<float> &posz, vector<float>& xx, vector<float>& yy, vector<float> & zz) {
    //cout<<comov_low<<" "<<comov_high<<" "<<comov_av<<" "<<comov_width<<endl;
    for ( int i=0 ;  i < posx.size(); i++){
       double radius= sqrt(pow(posx[i],2) + pow(posy[i],2) + pow(posz[i],2));
       //cout<<"radius= "<<radius<<endl;
       if(radius>= comov_low && radius <comov_high) {
	 xx.push_back(posx[i]); yy.push_back(posy[i]); zz.push_back(posz[i]); 
       }  
     }

    //posx.clear(); posy.clear(); posz.clear();
   }

  /* void calc_lightcone(vector<float>& xx, vector<float>& yy, vector<float>& zz) {
      map_a_particle_position_to_stacked_boxes();
      mapped_particles_in_snapshot_slice(xx, yy, zz); 
      }*/
  double get_comovinglow() {return comov_low;}
  double get_comovinghigh() {return comov_high;}
  double get_comovingav() {return comov_av;}
  double get_comovingwidth() {return comov_width;}  
  
  void print_pos(ofstream &outfile, vector<float> &x_lc, vector<float> &y_lc, vector<float> &z_lc) { 
    
    for ( int i=0 ; i < x_lc.size(); i++)
       outfile<<x_lc[i]<<" "<<y_lc[i]<<" "<<z_lc[i]<<endl;
    x_lc.clear(); y_lc.clear(); z_lc.clear();
    }

   void print_pos_binary(ofstream &outfile, vector<float> &x_lc, vector<float> &y_lc, vector<float> &z_lc) { 
     
     for ( int i=0 ; i < x_lc.size(); i++) {
      outfile.write(reinterpret_cast<const char*>(&x_lc[i]), sizeof (float));
      outfile.write(reinterpret_cast<const char*>(&y_lc[i]), sizeof (float));
      outfile.write(reinterpret_cast<const char*>(&z_lc[i]), sizeof (float));         }
      x_lc.clear(); y_lc.clear(); z_lc.clear();
    }

   void close_output_file(ofstream &outfile) {
    outfile.clear();  
    outfile.close();  
  }

 };
#endif
