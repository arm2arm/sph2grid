
//

#include <iostream>
#include <fstream>
#include <sstream> 
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include "data_readers.h"
#include "sphvolume.h"
using namespace std;

CTimer my_timer;
void usage(const char*execname, bool help=false)
{
if(help){
  cout<<"\n***********************************************************"<<endl;
  cout<<"A.Khalatyan"<<endl;
  cout<<"sph2grid v0.1 2009, Marseille"<<endl;
  cout<<"\n***********************************************************"<<endl;
  cout<<"Reads RHO and HSML blocks for the gas particles from the input file and\n"<<
    "generates density file on the regular grid defined by the input parameters.\n"<<
    "If user require Type>0 of the particles it will seek files with extension \"*_rho_type\".\n"<<
    "Example:\n"<<
    "If you have a snapshot name snap_200 and you require TYPE=1 then\n"<<
      "it will seek file snap_200_rho_1 with blocks \"H1\" for the smoothing lengths\n"<<
    "and \"RHO1\" for the densities.\n"<<endl;
  }

  cout<<"\nUsage:\n"<<execname<<" [FILEIN] [type] [XC YC ZC RC] [GRID]"<<endl;	
 
if(!help) cout<<"For more help:\n"<<execname<<" -help \n"<<endl;
 if(help) cout<<"***********************************************************\n"<<endl;

}


int main(int argc, char* argv[])
{
  if(argc!=8)
    {
      //Under linux directory separator given by "/"
      char* pch=strrchr(argv[0], '/');
      bool help=false;
      string execname;
      //Under windows directory separator given by "\"
      if(pch==NULL)
	{
	  pch=strrchr(argv[0], '\\');	  	    
	}
      if(pch!=NULL)
	execname=string(pch+1);
      else
	execname=string(argv[0]);
      
      if(argc==2)
	if((strstr(argv[1], "-h")!=NULL || strstr(argv[1], "/h")!=NULL ||strstr(argv[1], "?")!=NULL))
	  help=true;
      
      usage(execname.c_str(), help);
    }
  else
    {
      string fname=string(argv[1]);
      string foutname=fname;
      int isnap;
      if(!GetSnapName(foutname, isnap))
	{     
	  cout<<"Cannot extract file name from file: "<<fname<<endl;
	  return EXIT_FILE_NOT_FOUND;
	}
      foutname=foutname+"_grid";
      cout<<"In File: "<<fname<<endl;
      cout<<"Out File: "<<foutname<<endl;
      cout<<"Isnap: "<<isnap<<endl;
      /////////////////////////////////////////
      int type=atoi(argv[2]), GRID=atoi(argv[7]);
      float 
	XC=atoi(argv[3]),
	YC=atoi(argv[4]),
	ZC=atoi(argv[5]),
	RC=atoi(argv[6]);      
      /////////////////////////////////////////
      
      DoSph2Grid(fname, foutname,type,XC,YC,ZC,RC,GRID);
      /////////////////////////////////////////
      cout<<"done"<<endl;
    }

  return EXIT_SUCCESS;
}

