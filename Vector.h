// vectors: overloading operators example
#ifndef Vector_Header_included
#define Vector_Header_included

#include <iostream>
#include <vector>
using namespace std;

class CVector {
public:
  double x,y,z;
  CVector () {};
  
  CVector (double a, double b, double c, double par) {
     x = a*par;
     y = b*par;
     z = c*par;
  }

  CVector operator+ (CVector param) {
     CVector temp;
     temp.x = x + param.x;
     temp.y = y + param.y;
     temp.z=  z + param.z;
     return (temp);
   }
  CVector operator- (CVector param) {
    CVector temp;
    temp.x = x - param.x;
    temp.y = y - param.y;
    temp.z = z - param.z;
    return (temp);
  }
  CVector operator* (CVector param) {
    CVector temp;
    temp.x = x*param.x;
    temp.y = y*param.y;
    temp.z = z*param.z;
    return (temp);
  }
  friend class lightcone;
};





#endif


