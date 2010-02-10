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
#include <cstdlib>
#include <string>
#include "utils.h"


CRender::CRender():color_table(4),m_outfile("sph2grid.png") {
    // Just print the current CImg settings.
    //cimg_library::cimg::info();
    this->SetMinMaxRho(-7.0f,3.0f);

    std::string kisRender("SPH2GRID_DISPLAY");  
    char *pointer=NULL;
    show_image=false;
    if(pointer=getenv(kisRender.c_str()))show_image=true;

}

CRender::CRender(const CRender& orig) {
}

CRender::~CRender() {
}

void CRender::DumpAllLUT(void )
{
  for(int iLUT=0;iLUT<CImg<float>::get_MAX_NUM_LUT();iLUT++)
    {
      CImg<float> palette=CImg<float>::IDL_LUT256(iLUT);
      
      char buf[128];  
      sprintf(buf,"lut_%03d.png",iLUT);
      std::cout<<"Dumping LUT: "<<buf<<std::endl;
      palette.save(buf);
    }
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
    if(show_image)
      (image.normalize(0, 255)).display();



}

float Wsph(float u) {
    float retval = 0;
 
    if (0 <= u && u <= 1)
        retval = 1.0f - 3.0f / 2.0f * u * u + 3.0f / 4.0f * u * u * u;
    else
        if (1 < u && u <= 2)
        retval = 1.0f / 4.0f * (2.0f - u)*(2.0f - u)*(2.0f - u);
    else
        return 0;

    return retval;//(3.1456f*h*h*h);


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
    float h1=1.0f/hsml;
    for (int ip = 0; ip < np; ip++) {
        
        irange[0] = (int) (X[ip] - hsml);
        irange[1] = (int) (X[ip] + hsml);
        jrange[0] = (int) (Y[ip] - hsml);
        jrange[1] = (int) (Y[ip] + hsml);
        clamp2<int>(irange, 0, GRID);
        clamp2<int>(jrange, 0, GRID);
        
        //factor_normalize=1.0f/(4.0f/3.0f*M_PI*(hsml*hsml*hsml));
        for (j = jrange[0]; j < jrange[1]; j++) {
            dy=(j-Y[ip])*h1;
            for (i = irange[0]; i < irange[1]; i++)
            {
                dx=(i-X[ip])*h1;
                dist=(dx*dx+dy*dy);
                //if(dist>4.0f) continue;
                
                //std::cout<<"i="<<i<<" j="<<j<<" "<<irange[0]<<" "<<irange[1]
                //        <<" dx="<<dx<<" dy="<<dy<<" dist="<<dist<<std::endl;
                // std::cin.get();
                
                image(j, i, 0, 0) += kernel.W(std::sqrt(dist)); //Red
                //image(i, j, 0, 1) = image(i, j); //Green
                //image(i, j, 2) = image(i, j); //Blue
		range.getbound(image(j,i,0,0));

            }
            
        }
        if(ip % 1000==0){std::cout<<ip<<"     \r";std::cout.flush();}
    }
    range.print("Before scaling log:");
   
    // std::cout<<" "<<min_rho<<" "<<max_rho<<std::endl;
    range.Reset();
    for(int x=0;x<GRID;x++)
    for(int y=0;y<GRID;y++)
      {
        if (image(x, y, 0,0) > 0)
            image(x, y, 0,0) = std::log10(image(x, y, 0,0));
        else
            image(x, y, 0,0) =this->GetMinRho();
	image(x, y, 0,0) = clamp1<float>(image(x, y, 0,0),this->GetMinRho(), this->GetMaxRho());
       range.getbound(image(x,y,0,0));
    }

    range.print("pixel range:");

    
    image = image.normalize(0, 255); 
    CImg<unsigned char> palette=CImg<unsigned char>::IDL_LUT256(color_table);
    cimg_forXY(image, x, y) {
      unsigned int  ind=(unsigned int)image(x, y);
      const unsigned char col[]={palette(0,ind,0),palette(0,ind,1),palette(0,ind,2)};
      image(x,y,0,0)=col[0];
      image(x,y,0,1)=col[1];
      image(x,y,2)=col[2];

    }

    
    if(show_image)        
      (image).display();
    
    std::cout<<"Dumping PNG file"<<std::endl;
    image.save(m_outfile.c_str());


}

/// This routine uses adaptive size(HSML) for the partices
void CRender::DoRenderByAdaptivePoints(float *X, float *Y, float *Z,float *rho, float *hsml, int np, int GRID) {

    CImg<float> image(GRID, GRID, 1, 3,0);
    CRange range,range1;
   
    // loop over the all ip particles and find the image(i,j)
    // cells where it will contribute the mass.
    // the i and j are the image pixel coordinates;
    int i = 0, j = 0, irange[2], jrange[2];
    float  dx, dy, dist, factor_normalize=1.0;

    float *sqrtvec=new float[1000];
    for(int i=0;i<1000;i++)
      sqrtvec[i]=std::sqrt(i/1000.0*8.0);
    float *alpha=new float[256];
    for(int i=0;i<256;i++)
      {
	alpha[i] = 1.0;
	float w=1.0-i/256.0;
	for(int j=0;j<11;j++)
	  alpha[i] *= w;
	//alpha[i]*=0.5;
	
      }
      

    CEpanechikov<float> kernel(3);
    float h1=1.0f;
    for (int ip = 0; ip < np; ip++) {
        
        irange[0] = (int) (X[ip] - hsml[ip]);
        irange[1] = (int) (X[ip] + hsml[ip]);
        jrange[0] = (int) (Y[ip] - hsml[ip]);
        jrange[1] = (int) (Y[ip] + hsml[ip]);
        clamp2<int>(irange, 0, GRID);
        clamp2<int>(jrange, 0, GRID);
        h1=1.0f/hsml[ip];
	
        for (j = jrange[0]; j <= jrange[1]; j++) {
            dy=(j-Y[ip])*h1;
            for (i = irange[0]; i <= irange[1]; i++)
            {
                dx=(i-X[ip])*h1;
                dist=(dx*dx+dy*dy);
		//if(dist>8.0f) continue;
                
                //std::cout<<"i="<<i<<" j="<<j<<" "<<irange[0]<<" "<<irange[1]
                //        <<" dx="<<dx<<" dy="<<dy<<" dist="<<dist<<std::endl;
                // std::cin.get();
		//float val=kernel.W(std::sqrt(dist))*rho[ip];
		//	int ia= (int)(clamp1<float>(val,GetMinRho(), GetMaxRho())/(GetMaxRho()-GetMinRho())*255); 
		image(j, i, 0, 0) = image(j, i, 0, 0)+kernel.W(std::sqrt(dist))*rho[ip];
		/*alpha[ia]*///val; //Red

		//image(i, j, 0, 1) = image(i, j); //Green
                //image(i, j, 2) = image(i, j); //Blue
		//	range.getbound(image(j, i, 0, 0));
		//range1.getbound(dist);

            }
            
        }
        if(ip % 1000==0){std::cout<<ip<<"     \r";std::cout.flush();}
    }
    range1.print("DistRange:");
    range.print("Before scale the pixel range: ");
    range.Reset();
    cimg_forXY(image, x, y) {
        if (image(x, y, 0,0) > 0)
            image(x, y, 0,0) = std::log10(image(x, y, 0,0));
        else
            image(x, y, 0,0) =this->GetMinRho();
       image(x, y, 0,0) = clamp1<float>(image(x, y, 0,0),this->GetMinRho(), this->GetMaxRho());
       image(x, y, 0,1)=image(x, y, 0,0);
       image(x, y, 2)=image(x, y, 0,0);
       range.getbound(image(x,y,0,0));
    }

    range.print("pixel range:");

    
    image = image.normalize(0, 255);
   CImg<unsigned char> palette=CImg<unsigned char>::IDL_LUT256(color_table);
    cimg_forXY(image, x, y) {
      int ind=(int)image(x, y);
      const unsigned char col[]={palette(0,ind,0),palette(0,ind,1),palette(0,ind,2)};
      image(x,y,0,0)=col[0];
      image(x,y,0,1)=col[1];
      image(x,y,2)=col[2];

    }

    if(show_image)
        (image).display();
    
    std::cout<<"Dumping PNG file"<<std::endl;
    image.save(m_outfile.c_str());


}

/// This routine uses adaptive size(HSML) for the partices
void CRender::DoRenderByAdaptiveSortedPoints(float *X, float *Y,int *idx,float *rho, float *hsml, int np, int GRID) {

    CImg<float> image(GRID, GRID, 1, 3,0);
    CRange range,range1;
   
    // loop over the all ip particles and find the image(i,j)
    // cells where it will contribute the mass.
    // the i and j are the image pixel coordinates;
    unsigned int i = 0, j = 0, irange[2], jrange[2];
    float  dx, dy, dist;
    CEpanechikov<float> kernel(3);
    float h1=1.0f;
    int ip=0;
    unsigned int ia =0;
    float mir=std::pow(10.0f,min_rho),
      mar=std::pow(10.f,max_rho),
      scv=255.0f/(mar-mir)
      ; 
    
    float *alpha=new float[256];
    for(int i=0;i<256;i++)
      {
	alpha[i] = 1.0;
	float w=1.0f-i/256.0;
	for(int j=0;j<6;j++)
	  alpha[i] *= w;
	alpha[i] *= 1.0f;
	//alpha[i]=alpha[i]*(mar-mir)+mir;
      }
      

    for (int ipp = 0; ipp < np; ipp++) {
      ip=idx[ipp];
      irange[0] = (unsigned int) (X[ip] - hsml[ip]);
      irange[1] = (unsigned int) (X[ip] + hsml[ip]);
      jrange[0] = (unsigned int) (Y[ip] - hsml[ip]);
      jrange[1] = (unsigned int) (Y[ip] + hsml[ip]);
      clamp2<unsigned int>(irange, 0, GRID);
      clamp2<unsigned int>(jrange, 0, GRID);
      h1=1.0f/hsml[ip];
      
      for (j = jrange[0]; j <= jrange[1]; j++) {
	dy=(j-Y[ip])*h1;
	for (i = irange[0]; i <= irange[1]; i++)
	  {
	    dx=(i-X[ip])*h1;
	    dist=(dx*dx+dy*dy);
	    
	    float val=kernel.W(std::sqrt(dist))*rho[ip];
	    if(val>0)
	      {
		ia= (unsigned int)((val- mir)*scv); 
		image(j, i, 0, 0) = std::max(image(j, i, 0, 0),
					     /*alpha[ia]*/val);
		range.getbound(image(j, i, 0, 0));
	      }
	  }
	
        }
      range1.getbound(dist);
      if(ipp % 1000==0){std::cout<<ipp<<"     \r";std::cout.flush();}
    }
    range1.print("DistRange:");
    range.print("Before scale the pixel range: ");
    range.Reset();
    cimg_forXY(image, x, y) {
        if (image(x, y, 0,0) > 0)
            image(x, y, 0,0) = std::log10(image(x, y, 0,0));
        else
            image(x, y, 0,0) =this->GetMinRho();
       image(x, y, 0,0) = clamp1<float>(image(x, y, 0,0),this->GetMinRho(), this->GetMaxRho());
       image(x, y, 0,1)=image(x, y, 0,0);
       image(x, y, 2)=image(x, y, 0,0);
       range.getbound(image(x,y,0,0));
    }

    range.print("pixel range:");

    
    image = image.normalize(0, 255);
   CImg<unsigned char> palette=CImg<unsigned char>::IDL_LUT256(color_table);
    cimg_forXY(image, x, y) {
      int ind=(int)image(x, y);
      const unsigned char col[]={palette(0,ind,0),palette(0,ind,1),palette(0,ind,2)};
      image(x,y,0,0)=col[0];
      image(x,y,0,1)=col[1];
      image(x,y,2)=col[2];

    }

    if(show_image)
        (image).display();
    
    std::cout<<"Dumping PNG file"<<std::endl;
    image.save(m_outfile.c_str());


}
