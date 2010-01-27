/* 
 * File:   Render.cpp
 * Author: arm2arm
 * 
 * Created on January 26, 2010, 5:44 PM
 */

#include "Render.h"
#include <iostream>


/////////////////////// List of the CImg pluguins
#define cimg_verbosity 1
// The color tables:
#define cimg_plugin1 "plugins/LUT_ALL.rgb.h"

#include "CImg.h"
using namespace cimg_library;

#include <algorithm>
#include <cmath>
#include <istream>
#include "utils.h"


CRender::CRender() {
    // Just print the current CImg settings.
    cimg_library::cimg::info();
    this->SetMinMaxRho(-7.0f,3.0f);
}

CRender::CRender(const CRender& orig) {
}

CRender::~CRender() {
}

void CRender::DoRenderByGrid(float ***vol3d, int GRID, float zfac) {

    CImg<float> image(GRID, GRID, 1, 3);
    CRange range;
    for (int i = 0; i < GRID - 1; i++)
        for (int j = 0; j < GRID - 1; j++)
            for (int k = 0; k < GRID * zfac - 1; k++) {
                float pix = vol3d[i][j][k];
                if (pix < 1e-6)pix = 1.e-6;

                range.getbound(pix);
                image(i, j, 0, 0) = std::max(pix, image(i, j)); //Red
                image(i, j, 0, 1) = image(i, j); //Green
                image(i, j, 2) = image(i, j); //Blue

            }

    cimg_forXYC(image, x, y, v) {
        if (image(x, y, v) == 0) image(x, y, v) = 1e-6;
    }

    range.print("pixel range:");
    image = image.log10();
    image = image.normalize(image.mean(), image.max());
    (image.normalize(0, 255)).display();



}

template <class T>
void clamp(T *val, int num_elems, T min_value, T max_value) {
    for (int i = 0; i < num_elems; i++) {
        if (val[i] < min_value)val = min_value;
        else
            if (val[i] > max_value)val = max_value;
    }
}

template <class T>
void clamp2(T(&val)[2], T min_value, T max_value) {
    if (val[0] < min_value)val[0] = min_value;
    else
        if (val[0] > max_value)val[0] = max_value;

    if (val[1] < min_value)val[1] = min_value;
    else
        if (val[1] > max_value)val[1] = max_value;
}
template <class T>
T clamp1(T val, T min_value, T max_value) {
    if (val < min_value)return  min_value;
    else
        if (val > max_value)return  max_value;
    return val;
}
//Smoothing kernel

float kernel(float u) {
    //standard cubic spline
    if (u < 1.0)
        return (1. - 1.5 * u * u + 0.75 * u * u * u);
    else if (u < 2.0)
        return 0.25 * (2. - u)*(2. - u)*(2. - u);
    else
        return 0.0f;
}

//////////////////////////////////////
//X , Y and Z must be in a grid units.

void CRender::DoRenderByPoints(float *X, float *Y, float *Z, float hsml, int np, int GRID) {

    CImg<float> image(GRID, GRID, 1, 3,0);
    CRange range;

    // loop over the all ip particles and find the image(i,j)
    // cells where it will contribute the mass.
    // the i and j are the image pixel coordinates;
    int i = 0, j = 0, irange[2], jrange[2];
    float  dx, dy, dist, factor_normalize=1.0;
    
    CEpanechikov<float> kernel(3);
    for (int ip = 0; ip < np; ip++) {
        irange[0] = (int) (X[ip] - hsml);
        irange[1] = (int) (X[ip] + hsml);
        jrange[0] = (int) (Y[ip] - hsml);
        jrange[1] = (int) (Y[ip] + hsml);
        clamp2<int>(irange, 0, GRID);
        clamp2<int>(jrange, 0, GRID);
        //factor_normalize=1.0f/(4.0f/3.0f*M_PI*(hsml*hsml*hsml));
        for (j = jrange[0]; j < jrange[1]; j++) {
            dy=(j-Y[ip])/hsml;
            for (i = irange[0]; i < irange[1]; i++)
            {
                dx=(i-X[ip])/hsml;
                dist=std::sqrt(dx*dx+dy*dy);
            /*    std::cout<<"i="<<i<<" j="<<j<<" "<<irange[0]<<" "<<irange[1]
                        <<" dx="<<dx<<" dy="<<dy<<" dist="<<dist<<std::endl;
                 std::cin.get();
                */
                image(i, j, 0, 0) += kernel.W(dist); //Red
                //image(i, j, 0, 1) = image(i, j); //Green
                //image(i, j, 2) = image(i, j); //Blue


            }
            
        }
        if(ip % 100==0){std::cout<<ip<<"     \r";std::cout.flush();}
    }

    cimg_forXY(image, x, y) {
        if (image(x, y, 0,0) > 0)
            image(x, y, 0,0) = std::log10(image(x, y, 0,0));
        else
            image(x, y, 0,0) =this->GetMinRho();
       image(x, y, 0,0) = clamp1<float>(image(x, y, 0,0),this->GetMinRho(), this->GetMaxRho());
       range.getbound(image(x,y,0,0));
    }

    range.print("pixel range:");

    
    image = image.normalize(0, 255);
   CImg<float> palette=CImg<float>::IDL_LUT256(13);
    cimg_forXY(image, x, y) {
      int ind=image(x, y);
      const unsigned char col[]={palette(0,ind,0),palette(0,ind,1),palette(0,ind,2)};
      image(x,y,0,0)=col[0];
      image(x,y,0,1)=col[1];
      image(x,y,2)=col[2];

    }

    if(GRID <800)
        (image).display();
    
    std::cout<<"Dumping PNG file"<<std::endl;
    image.save("sph2grid.png");


}

