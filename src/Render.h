/* 
 * File:   Render.h
 * Author: arm2arm
 *
 * Created on January 26, 2010, 5:44 PM
 */

#ifndef _RENDER_H
#define	_RENDER_H
#include <string>
#include <iostream>
class CRender {
 public:
  CRender();
  CRender(const CRender& orig);
  virtual ~CRender();
  void DoRenderByGrid(float ***vol3d, int GRID=128, float zfac=1.0);
  void DoRenderByPoints(float *X, float *Y,float *Z, float hsml, int np, int GRID=128);
  void DoRenderByAdaptivePoints(float *X, float *Y, float *Z,float *rho, float *hsml, int np, int GRID);
  void DoRenderByAdaptiveSortedPoints(float *X, float *Y,int *idx,float *rho, float *hsml, int np, int GRID);
  void DoSPHVolume(float *X, float *Y, int *idx,float *rho, float *hsml, int np, int GRID);

  float GetMinRho(void){return min_rho;};
  float GetMaxRho(void){return max_rho;};
  void SetMinMaxRho(float _min_rho, float _max_rho){
    min_rho=_min_rho;max_rho=_max_rho;
  };
  void DumpAllLUT(void );
  void   SetColorTable(int i){color_table=i-1;
  std::cout<<"Loading color table="<<color_table<<std::endl;};
  void SetOut(std::string s){m_outfile=s;}
 private:
  float  Wsph(float, float);
  float min_rho, max_rho;
  int color_table;
  bool show_image;
  std::string m_outfile;
};

#endif	/* _RENDER_H */

