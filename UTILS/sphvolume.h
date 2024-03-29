#ifndef _SPHVOLUME_
#define _SPHVOLUME_
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace std;
/****************/
void AllocateVol3D();
float Wsph(float rr, float h);/*kerlen for the SPH*/
void DeallocateVol3D();
void AllocateVol3D();
int WriteSPHVolume(string fname);
void  DoSPHVolume();/* use Rho and HSML from RAW file and assign to the 3D volume grid*/
void  SmoothSPHVolume(int sm ); /*BoxCar smoothing of the 3D volume, similar thing as IDL routine smooth*/
/****************/

#endif
