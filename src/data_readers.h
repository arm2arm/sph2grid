#ifndef _DATA_READERS_ 
#define _DATA_READERS_
//#pragma warning(disable:4996)

#include <functional>
#include <algorithm>
#include <string>
#include <utility>
#include <map>
#include <stack>
#include <vector>
#include <list>
#include <iostream>
#include <sstream> 
#include <stdlib.h>
#ifdef WIN32
#include <conio.h>
#endif
#include <stdio.h>
#include <ctime>
#include <ctype.h>
#include <list>
#include <cmath>
// Needed for the ifstream class
#include <fstream>
// Needed for the setw stream manipulator
#include <iomanip>
#ifdef MYTREE
#include "tree.h"
#endif
using namespace std;

/////////////////////////////////////////////////
#define  int4bytes int
#define EXIT_FILE_NOT_FOUND 13
/////////////////////////////////////////////////
float mydrand48(void);
bool FileExists(char* filename);
bool FileExists(const char* & filename);
bool FileExists(string filename);
bool StringToInt(const string &s, int &i);
bool GetSnapPath( string snap, string &path);
bool GetSnapName( string &snap, int &isnap);
class CTimer
	{
	public:
		void start(){time (&m_start);};
		void stop(bool print_flag=false){
			time (&m_end);
			m_dif = difftime (m_end,m_start);
			if(print_flag)
				print();
			}
		void print(void ){cout<<m_dif<<" secs"<<endl;};
		time_t m_start;
		time_t m_end;
		double m_dif;
	};
//////////////////////////////////////////////////////////
class CRegion
	{
	public:
		bool PointInRect(float *pt)
			{
			if(PERIODIC_FLAG)
				{
				pt[0]-=tx;
				pt[1]-=ty;
				pt[2]-=tz;

				pt[0]=periodic( pt[0], xmax,xmin);
				pt[1]=periodic( pt[1], xmax,xmin);
				pt[2]=periodic( pt[2], xmax,xmin);
	 
				}
			if( pt[0]>(m_cent[0]-m_R) && pt[0]<(m_cent[0]+m_R)&&
				pt[1]>(m_cent[1]-m_R) && pt[1]<(m_cent[1]+m_R)&&
				pt[2]>(m_cent[2]-m_R) && pt[2]<(m_cent[2]+m_R))
				return true;
			return false;				
			};
		bool PointInSphere(float *pos)
			{
			 
			 float dist2=(pos[0]-m_cent[0])*
				 (pos[0]-m_cent[0])+
				 (pos[1]-m_cent[1])*
				 (pos[1]-m_cent[1])+
				 (pos[2]-m_cent[2])*
				 (pos[2]-m_cent[2]);
			 if((m_R*m_R)<dist2)
				 return true;
			 return false;
			
			};
		void SetPeriodic(float x, float y, float z,
			float xxmin, float xxmax)
			{
			 PERIODIC_FLAG=true;
			 xmin=xxmin;xmax=xxmax;
			 tx=x;ty=y;tz=z;
			};
		float periodic(float xp, float xmax, float xmin)
			{

			xp+=(xmax-xmin)/2.0f;
			if(xp>xmax)
				{
				xp=xp-(xmax-xmin);
				}

			if(xp<xmin) 
				{
				xp=xp+(xmax-xmin);
				}

			return xp;
			}
		CRegion(){};
		CRegion(float x,float y,float z, float R):tx(0.0f),ty(0.0f),tz(0.0f){
			m_R=R;
			xmin=0.0f; xmax=0.0f;
			m_cent[0]=x;
			m_cent[1]=y;
			m_cent[2]=z;
			PERIODIC_FLAG=false;
			};
		float m_R;
		float tx, ty, tz;
		float m_cent[3];
		float xmin, xmax;
		bool PERIODIC_FLAG;
	       
	};
//////////////////////////////////////////////////////////
struct strParticleData
	{
	float Pos[3];
	};
struct strSPHData
	{
	  float Rho;
	  float Hsml;
	  float Temp;
	  float EnDt;
	  float Type;
	};
struct strSPHParticle
	{
	 float Pos[3];
	 strSPHData sph;
	 strSPHParticle(){};
	 strSPHParticle(float *x,char inType)
		{
		memcpy(&Pos[0], x, sizeof(float)*3);
		sph.Type=inType;
		};
	strSPHParticle(float *x, float rho=0,float hsml=0, float temp=0)
		{
		memcpy(&Pos[0], x, sizeof(float)*3);
		sph.Rho=rho;sph.Hsml=hsml;sph.Temp=temp;
		sph.Type=0;
		};
	};
class strSPHParticle_cmp
	: std::binary_function<const strSPHParticle&, const strSPHParticle&, bool>
	{
	public:
	enum DATASORTMODE
		{
		AscRho,
		DescRho,
		AscTemp,
		DescTemp
		};

	strSPHParticle_cmp(DATASORTMODE mode=AscRho):m_mode(mode){};
	bool operator() (const strSPHParticle& lhs, const strSPHParticle& rhs ) const
		{
		if(m_mode==DescRho)
			return lhs.sph.Rho<rhs.sph.Rho;
		if(m_mode==AscRho)
			return lhs.sph.Rho>rhs.sph.Rho;
		if(m_mode==DescTemp)
			return lhs.sph.Temp<rhs.sph.Temp;
		if(m_mode==AscTemp)
			return lhs.sph.Temp>rhs.sph.Temp;
		return lhs.Pos[2]<rhs.Pos[2];
		}
	private:
		DATASORTMODE m_mode;
	};
class CAsciiReader
	{
	public:
		CAsciiReader(){};
		~CAsciiReader(){};
		CAsciiReader(string IN_FILE, unsigned int nfields=1):m_filename(IN_FILE),m_nfields(nfields){
			m_good=ReadAsciiFile();
			};
		typedef vector<float> dynvector;
		map<unsigned int,dynvector> m_data;
		bool ReadAsciiFile()
			{
			if(!FileExists(m_filename))
				{
				cout<<"Cannot find File..."<<endl;
				cout<<"Failed to open: "<<m_filename<<endl;
				return false;
				};
			ifstream ifile(m_filename.c_str());
			m_data.clear();
			while(!ifile.eof())
				{
				float a;
				int snap;
				ifile>>snap;
				for(unsigned int i=0;i<m_nfields;i++)
					{
					ifile>>a;
					m_data[snap].push_back(a);
					}
				}
			ifile.close();
			return true;
			};
		bool IsGood(){return m_good;};
	protected:
		bool m_good;
		string m_filename;
		unsigned int m_nfields;
	};

class CData
	{
	public:
		CData(){ID=NULL;P=NULL;};
	virtual	~CData()
			{
			//cout<<"CData::~CData freeding memory";
			 FreeMemory();
			//cout<<"..Ok"<<endl;
			};
		void FreeMemory()
			{
			if(P!=NULL)
			 delete [] P;
			if(ID!=NULL)
			 delete [] ID;
			 m_isgood=false;
			 m_data.clear();
			};
		bool good(){return m_isgood;};
		int *ID;
		strParticleData *P;
		strSPHData *sph;
		vector<strSPHParticle> m_data;
		void SortParticlesByRho(bool mode=true)
			{
			if(mode)
			sort(m_data.begin(), m_data.end(), 
				strSPHParticle_cmp(strSPHParticle_cmp::DescRho));
			else
			sort(m_data.begin(), m_data.end(), 
				strSPHParticle_cmp(strSPHParticle_cmp::AscRho));
			};
		void SortParticlesByTemp(bool mode=true)
			{
			if(mode)
			sort(m_data.begin(), m_data.end(), 
			strSPHParticle_cmp(strSPHParticle_cmp::DescTemp));
			else
			sort(m_data.begin(), m_data.end(), 
				strSPHParticle_cmp(strSPHParticle_cmp::AscTemp));
			};
		void ScaleVect(float *v, unsigned int nelem, float xmin, float xmax, float vBox)
			{
				for(unsigned int i=0;i<nelem;i++)
					{
					v[i]=(v[i]-xmin)/(xmax-xmin);
					v[i]=v[i]*2.0f*vBox-vBox;
					}
			};
		void ScaleVectNoMove(float *v, unsigned int nelem, float xmin, float xmax, float vBox)
			{
				for(unsigned int i=0;i<nelem;i++)
					{
					v[i]=(v[i]-xmin)/(xmax-xmin);
					v[i]=v[i]*vBox;
					if(v[i]<=0)
						v[i]=1;
					else
						if(v[i]>=vBox)v[i]=vBox;
					}

			};
		void PutInCOM(strParticleData *P, unsigned int n, float  *COM_FROM_FILE=NULL)
			{
			
			unsigned int i;
			if(COM_FROM_FILE==NULL)
			{
			memset(m_COM, 0, sizeof(m_COM));
			for(i=0;i<n;i++)
				{
				m_COM[0]+=P[i].Pos[0];
				m_COM[1]+=P[i].Pos[1];
				m_COM[2]+=P[i].Pos[2];
				}
			m_COM[0]/=float(n);
			m_COM[1]/=float(n);
			m_COM[2]/=float(n);
			}else
				{
				memcpy(m_COM, COM_FROM_FILE, sizeof(float)*3);				
				}
			for(i=0;i<n;i++)
				{
				P[i].Pos[0]-=m_COM[0];
				P[i].Pos[1]-=m_COM[1];
				P[i].Pos[2]-=m_COM[2];
				}

			};
		void PutInCOM(float  *COM_FROM_FILE=NULL)
			{
			unsigned int n=m_data.size();
			unsigned int i;
			if(COM_FROM_FILE==NULL)
			{
			memset(m_COM, 0, sizeof(m_COM));
			for(i=0;i<n;i++)
				{
				m_COM[0]+=m_data[i].Pos[0];
				m_COM[1]+=m_data[i].Pos[1];
				m_COM[2]+=m_data[i].Pos[2];
				}
			m_COM[0]/=float(n);
			m_COM[1]/=float(n);
			m_COM[2]/=float(n);
			}else
				{
				memcpy(&m_COM[0], &COM_FROM_FILE[0], sizeof(float)*3);				
				}
			for(i=0;i<n;i++)
				{
				m_data[i].Pos[0]-=m_COM[0];
				m_data[i].Pos[1]-=m_COM[1];
				m_data[i].Pos[2]-=m_COM[2];
				}

			};
		void AutoCenterByRho()
			{
			
				unsigned int i, Ntrace;
				memset(m_COM, 0, sizeof(m_COM));
				SortParticlesByRho(false);
				Ntrace=(unsigned int)(0.1*m_data.size());
				for(i=0;i<Ntrace;i++)
					{
					m_COM[0]+=m_data[i].Pos[0];
					m_COM[1]+=m_data[i].Pos[1];
					m_COM[2]+=m_data[i].Pos[2];
					}
				m_COM[0]/=float(Ntrace);
				m_COM[1]/=float(Ntrace);
				m_COM[2]/=float(Ntrace);
				for(i=0;i<(unsigned int)m_data.size();i++)
					{
					m_data[i].Pos[0]-=m_COM[0];
					m_data[i].Pos[1]-=m_COM[1];
					m_data[i].Pos[2]-=m_COM[2];
					}
				cout<<"\n----------------\n"<<
					"Center of mass: "<<m_COM[0]<<", "
					<<m_COM[1]<<", "<<m_COM[2]<<
					"\n----------------\n"<<endl;
			};

		float m_COM[3];
		void SetDataName(string name){m_dataname=name;};
		string GetDataName(void){return m_dataname;};
		vector<strSPHParticle> m_bhdata;
	protected:
		bool m_isgood;
		string m_dataname;
		
	};

class CReaders :public CData
	{
	public:
		CReaders():swp_flag(false),find_flag(false)
			{};
		virtual ~CReaders() {};
		virtual bool ReadData(string file){return true;};
		int GetBlk(ifstream *file_to_read, int *blk)
			{

			file_to_read->rdbuf()->sgetn(reinterpret_cast<char *>(blk), (std::streamsize)sizeof(int));
			m_isgood=file_to_read->bad();
			return *blk;
			};
		void GetBlkName(ifstream *file_to_read, char *name)
			{

			file_to_read->rdbuf()->sgetn(name,sizeof(char)*4);
			GetBlk(file_to_read, &blk);
			};

		void swap_Nbyte(char *data,int n,int m)
			{
			int i,j;
			char old_data[16];

			if(swp_flag)
				{
				for(j=0;j<n;j++)
					{
					memcpy(&old_data[0],&data[j*m],m);
					for(i=0;i<m;i++)
						{
						data[j*m+i]=old_data[m-i-1];
						}
					}
				}
			};
		size_t my_fread(void *ptr, size_t size, size_t nmemb, ifstream *file_to_read)
			{
			size_t nread;

			if((nread = file_to_read->rdbuf()->sgetn((char*)ptr,size*nmemb)) != size*nmemb)
				{
				if(!find_flag)
					{
					cerr<<"I/O error (fread) "<<nread<<endl;
					exit(3);
					}else
					{
					cerr<<"Cannot find block: ";
					return 0;
						}
				}
			return nread;
			};

		bool swp_flag;
		bool find_flag;
		int blk;
		string m_filename;
		ifstream m_file;
	};


//////////////////////////////////////////////////////////
class CGadget : public  CReaders
	{
	public:
		struct io_header
			{
			int npart[6];
			double mass[6];
			double time;
			double redshift;
			int flag_sfr;
			int flag_feedback;
			int npartTotal[6];
			int flag_cooling;
			int num_files;
			double BoxSize;
			double Omega0;
			double OmegaLambda;
			double HubbleParam;
			int flag_multiphase;
			int flag_stellarage;
			int flag_sfrhistogram;
			int flag_metals;
			int flag_decouple;
			int flag_effmodel;
			char fill[72];		/* fills to 256 Bytes */
			}myhead;
		CGadget(string file, bool readnow=true):m_NBH(0)
			{
			m_filename=file;
			cout<<"Opening: "<< m_filename<<endl;
			if(!FileExists(m_filename)){cout<<"Can not find file: \n"<< m_filename<<endl;
				m_isgood=false;
			//exit(EXIT_FILE_NOT_FOUND);
				};
			m_file.open(m_filename.c_str(),  ios::in|ios::binary);

			m_isgood=GetFileFormat();
			if(readnow)
				{
				m_isgood=ReadData(m_filename);
				}
			GetHeader();
			rhoname.push_back("RHO ");
			rhoname.push_back("RHO1");
			rhoname.push_back("RHO2");
			rhoname.push_back("RHO3");
			rhoname.push_back("RHO4");
			rhoname.push_back("RHO5");

			hsmlname.push_back("HSML");
			hsmlname.push_back("H1  ");
			hsmlname.push_back("H2  ");
			hsmlname.push_back("H3  ");
			hsmlname.push_back("H4  ");
			hsmlname.push_back("HS5 ");
			
			uname.push_back("U   ");
			uname.push_back("RHO1 ");
			uname.push_back("RHO2 ");
			uname.push_back("RHO3 ");
			uname.push_back("AGE ");
			uname.push_back("MBH ");

			};
		int  unit_conversion(void);
		~CGadget(){
			indexinreg.clear();
			//cout<<"Deleting data in CGadget::~CGadget"<<endl;
			};
		bool ReadData(string file);
		void GetHeader(ifstream *file_to_read);
		void GetHeader(void){GetHeader(&m_file);};
		void GetHead(ifstream *fd, io_header &head);
		int find_block(ifstream *fd,const char *label);
		bool good(){return m_isgood;};	
		void SeekToType(ifstream *file_to_read,int type, int onesize, bool mass_flag=false);
		bool GetFileFormat();
		bool GetFileFormat(ifstream &filein);
		bool GetGas(CRegion reg, bool flag_putin_COM=false);
		bool GetMYGas(CRegion reg, int Type);
		bool GetStars(CRegion reg, bool flag_putin_COM=false);
		bool GetBH(CRegion reg);
		bool CheckRhoFile(CRegion reg);
		unsigned int read_block(float *&pV, char *name, int t);
		unsigned int read_blockv3(float *&pV, char *name, int t);
		void WriteRhoFile(string rhofilename, int type=4);
		void WriteOneBlock(ostream &file,string blname, char* pData, unsigned int datasize);
		bool GetSPHParticles(int type, CRegion reg, bool flag_putin_COM);
		bool m_isgood;
		string m_name;
		vector<unsigned int> indexinreg;
	public:
		vector<string> rhoname, hsmlname,uname;
	protected:
#ifdef MYTREE
	stack<CStone> m_matter;
#endif
	float Rhosc, Rhomax,Rhomin;
	public:
	int m_NBH;
	};
/////////////////// Typedefs for Galaxy sorting //////
class CGalaxy;

typedef pair <vector<CGalaxy*>::iterator, float> Map_Int_Flt_Pair;

typedef pair <int, int> Map_Int_Int_Pair;

////////////////////CGalaxy class ////////
class CGalaxy{
public:
	CGalaxy(int idgrp=0):grpID(idgrp){
		memset(x,0,sizeof(x));
		//		memset(Rxyz,0,sizeof(Rxyz));
		memset(npart,0,sizeof(npart));
		Ntotal=0;R=0;ID=0;m_Rxx=0;		
		};

	~CGalaxy(){};
	void SetIds(int nids){id.resize(nids,0);};
	void SortIds()
		{
		if(id_rxx.size()>0)
			sort(id_rxx.begin(), id_rxx.end());
		if(id.size()>0)
			{
			sort(id.begin(), id.end());
			int nold=id.size();
			id.erase (unique(id.begin(),
				id.end()),
				id.end());
			if(nold!=id.size())
				{
				cout << "Corrupted ID list for IDGrp="<<this->ID<<" "<<this->grpID<<" found "<<
					nold-id.size()<<" duplicates out of "<<nold<<" in ID list: "<< endl;
				//exit(0);
				}


			}
		};
	/////////////////////////////
	inline bool operator == (const CGalaxy &b) const
		{
		return (b.grpID==grpID);
		}
	/////////////////////////////
	//vector<CIntersect> child_list;
	vector<Map_Int_Flt_Pair> child_list;
	vector<int> id;
	vector<int> id_rxx;
	float x[3];
	float R;
	float m_Rxx;
	int ID;
	int grpID;
	int npart[3];
	//int Rxyz[3];
	int Ntotal;


	};

//////////////////: AFOF CLASS /////////////////
class CFOFCatalog : public CReaders
	{
	public:
		CFOFCatalog(){IDs=NULL;m_gadobj=NULL;};
		CFOFCatalog(string file);
		~CFOFCatalog(){
			vector<CGalaxy*>::iterator ithalos=halos.begin();
			while(ithalos!=halos.end())
				{
					//cout<<"~CFOFCatalog: deleting"<<endl;
					delete  (*ithalos);
					//ithalos=halos.erase(ithalos);
					++ithalos;
				};
			halos.clear();

			};
		bool ReadCatalog(string file);
		bool ReadCatalogIDs(string file);
		bool ReadCatalogIDsBin(string file);
		bool InitIDs( CGadget *gad)
			{
			Redshift=(float)gad->myhead.redshift;
			m_gadobj=gad;
			return InitIDs( gad->ID);
			}
		bool InitIDs( int *ID)
			{
			cout<<"Init IDS ...";
			bool ret=false;
			string file=m_afofcat, filebin;
			file.resize(m_afofcat.size()-4);
			IDs=ID;
			filebin=file;
			filebin+=string(".bin");
			cout<<".RB.";
			ret=ReadCatalogIDsBin(filebin);
			if(!ret)
				{
				cout<<".RA.";
				ret=ReadCatalogIDs(file);
				}
			cout<<"ok"<<endl;
			return ret;
			}
		void FillIDs();

		vector<CGalaxy*> halos;
		vector<unsigned int> grpIDs;
		int *IDs;
		CGadget* m_gadobj;
		float Redshift;

		vector<int> indexVector;// used for sorting 
		vector<int>::iterator it;// used for accessing to sorted elements
		string m_afofcat;

	};
////////////////////////////////////////////////
bool ReadAFOFHaloes(void);
bool ReadASCIIParticles(); // Read ASCII files
bool ReadGadgetParticles(bool first_flag=true, float tx=0, float ty=0, float tz=0); // Read Gadget  files
bool ReadARTParticles(); //Read ART Files
bool ReadMLAPMParticles();// MLAPM
			

#endif


