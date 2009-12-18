#include "data_readers.h"
#include "utils.h"
#define DEBUG_ME false
#define RadFrac  0.6
/////////////////////////////
CTimer timer;//timer for tests
//////////////////////
// File Utils ////////
/////////////////////////////
// Checks if a file exists //
/////////////////////////////
bool FileExists(char* filename) {

	ifstream file;
	file.open(filename);

	// Check if the file exists
	if (file.is_open() == true) {
		file.close();
		return true;
		}

	return false;

	}


bool FileExists(const char* & filename) {

	ifstream file;
	file.open(filename);

	// Check if the file exists
	if (file.is_open() == true) {
		file.close();
		return true;
		}

	return false;

	}


bool FileExists(string filename) {

	ifstream file;
	file.open(filename.c_str());

	// Check if the file exists
	if (file.is_open() == true) {
		file.close();
		return true;
		}

	return false;

	}


//////////////////////
bool StringToInt(const string &s, int &i)
	{
	istringstream myStream(s);

	if (myStream>>i)
		return true;
	else
		return false;
	}
bool GetSnapPath( string snap, string &path)
	{
	basic_string <char>::size_type indexCh;
	static const basic_string <char>::size_type npos = 0;

	indexCh=snap.find_last_of('/');
	if((indexCh)!=npos)
		path.assign(snap,0, int(indexCh));
	else
		path.assign("./");
	return true;
	}
bool GetSnapName( string &snapio, int &isnap)
	{
	string  dig_snap;
	string snap=snapio;
	basic_string <char>::size_type indexCh;
	static const basic_string <char>::size_type npos = 0;

	indexCh=snap.find_last_of('/');
	if((indexCh)!=npos)
		snap.assign(snap,int(indexCh+1), int(snap.size()-1));
	snapio=snap;
	if(strstr(snap.c_str(), "_")!=NULL)
	while(snap.find_first_of("_")!=0)
		{
		indexCh=snap.find_first_of("_");
		dig_snap=dig_snap.assign(snap,int(indexCh+1), int(3));
		if (StringToInt(dig_snap, isnap))
			{
			snap=dig_snap;
			return true;
			}
		else
			snap.assign(snap,int(indexCh+1), int(snap.size()-1));

		}
	return false;
	}
//////////////////////


size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
	{
	size_t nread;

	if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
		{
		printf("I/O error (fread) %d !\n",(unsigned int )nread);
		exit(3);
		}
	return nread;
	}


void CGadget::SeekToType(ifstream *file_to_read,int type, int onesize, bool  mass_flag)
	{
	long nskip=0;
	for(int ip=0;ip<type;ip++)
		{
		if(!mass_flag)
			nskip+=(this->myhead.npart[ip]);
		else
			{
			bool flag=(myhead.npart[ip]*myhead.mass[ip])==0.0?true:false;
			nskip+=(this->myhead.npart[ip]*flag);
			}
		}
	file_to_read->seekg((nskip)*onesize, ios_base::cur);
	}

int CGadget::find_block(ifstream *fd,const char *label)
	{
	int4bytes blocksize=0, blksize=0, ret=0;
	char blocklabel[5]={"    "};
	find_flag=true;
	fd->seekg(ios::beg);
	printf("Finding: %s\n",label);
	while(!fd->eof() && (blocksize == 0))//&& strcmp(blocklabel, "Z   ")!=0 )
		{
		//      cout<<"SKIP1"<<endl;
		GetBlk(fd, &blksize);
		//      cout<<"SKIP2"<<endl;
		//				if(blksize == 134217728)
			{				
			swap_Nbyte((char*)&blksize,1,4);
			}
			if(blksize != 8)
				{
				if(ret>0)
					{
					printf("incorrect format (blksize=%d)!\n",blksize);
					exit(1);
					}else break;
				}
			else
				{

				ret=my_fread(blocklabel, 4*sizeof(char), 1, fd);
				if(ret>0)
					ret=my_fread(&blocksize, sizeof(int4bytes), 1, fd);
				else
					return false;
				swap_Nbyte((char*)&blocksize,1,4);
				if( DEBUG_ME)
					{
					if(ret>0)
						printf("Found Block <%s> with %d bytes\n",blocklabel,blocksize);
					else						
						printf(" <%s> \n",label);


					}

				GetBlk(fd, &blksize);
				if(strcmp(label,blocklabel)!=0)
					{ 	
					fd->seekg(blocksize,ios_base::cur);
					blocksize=0;
					}
				}
		}
	find_flag=false;
	return(blocksize-8);
	}

bool CGadget::GetGas(CRegion reg,bool flag_putin_COM)
	{
	unsigned int ngas=myhead.npart[0], ninreg=0;
	unsigned int sizeall, size, ip;
	float *pF;
	CRange range, range1;
	sizeall=find_block(&m_file, "POS ");
	cout<<"POS"<<endl;
	GetBlk(&m_file, &blk);
	size=ngas*sizeof(float);
	P=new strParticleData [ngas];
	my_fread(P, size, 3, &m_file);
	swap_Nbyte((char*)P,ngas*3,4);
	if(flag_putin_COM)
		PutInCOM(P, ngas);
range.Reset();
for(ip=0;ip<ngas;ip++)
{
	range.getbound(P[ip].Pos[0]);
	range.getbound(P[ip].Pos[1]);
	range.getbound(P[ip].Pos[2]);
}
range.print("Pos Bounds: ");

	for(ip=0;ip<ngas;ip++)
		{
		if(reg.PointInRect(P[ip].Pos))
			{	
			indexinreg.push_back(ip);
			m_data.push_back(
				strSPHParticle(&P[ip].Pos[0]));
			}
		}
	ninreg=m_data.size();	

	sizeall=find_block(&m_file, "RHO ");
	cout<<"RHO"<<endl;
	GetBlk(&m_file, &blk);
	my_fread(P, size, 1, &m_file);
	pF=(float *)P;
	swap_Nbyte((char*)pF,ngas,4);
	for(ip=0;ip<ninreg;ip++)
		{
		 m_data[ip].sph.Rho=	pF[indexinreg[ip]];
		}
	sizeall=find_block(&m_file, "U   ");
	cout<<"U"<<endl;
	GetBlk(&m_file, &blk);
	my_fread(P, size, 1, &m_file);
	pF=(float *)P;
	swap_Nbyte((char*)pF,ngas,4);
	for(ip=0;ip<ninreg;ip++)
		{
		 m_data[ip].sph.Temp=	pF[indexinreg[ip]];
		}
	
	sizeall=find_block(&m_file, "HSML");
	cout<<"HSML"<<endl;
	GetBlk(&m_file, &blk);
	my_fread(P, size, 1, &m_file);
	pF=(float *)P;
	swap_Nbyte((char*)pF,ngas,4);
	for(ip=0;ip<ninreg;ip++)
		{
		 m_data[ip].sph.Hsml=	pF[indexinreg[ip]];
		}

	sizeall=find_block(&m_file, "ENDT");
	cout<<"ENDT"<<endl;
	GetBlk(&m_file, &blk);
	my_fread(P, size, 1, &m_file);
	pF=(float *)P;
	swap_Nbyte((char*)pF,ngas,4);
	range.Reset();
	range1.Reset();
	for(ip=0;ip<ninreg;ip++)
		{
		  m_data[ip].sph.EnDt=log10(abs(pF[indexinreg[ip]])+1)*
		    log10(m_data[ip].sph.Rho+1);
		 range.getbound(m_data[ip].sph.EnDt);
		}
	range.print("Entropy range:");
	unit_conversion();
	for(ip=0;ip<ninreg;ip++)
	  {
	    m_data[ip].sph.EnDt=(m_data[ip].sph.EnDt-range.Min)/
	      (range.Max-range.Min);
	    range1.getbound(m_data[ip].sph.EnDt);	    
	  }
	range1.print("Entropy LOG range:");

	delete P;
	P=NULL;
	indexinreg.clear();
	return true;
	};

// Check rhoi file if exist then read else make rho
// the rho file is inside the file snapfile+"_rho"
bool CGadget::CheckRhoFile(CRegion reg)
	{
	string rhofilename=m_filename+"_rho";
	if(!FileExists(rhofilename)){cout<<"Can not find Rho file to get density: \n"<< rhofilename<<endl;
	return false;			
	};
	
	return true;
	}
bool CGadget::GetStars(CRegion reg, bool flag_putin_COM)
	{
	return GetSPHParticles(4, reg,flag_putin_COM);	
	}

bool CGadget::GetMYGas(CRegion reg, int Type)
	{
	  return GetSPHParticles(Type, reg,false);//flag_putin_COM);	
	}
bool CGadget::GetBH(CRegion reg)
	{
	if(this->myhead.npart[5]>10)
		;//GetSPHParticles(5, reg,false);	
	else
		{
		unsigned int nbh=myhead.npart[5];	
		unsigned int sizeall, size, ip;
		float *pM=NULL;
		strParticleData *pVel=NULL;
		cout<<"POS"<<endl;
		sizeall=find_block(&m_file, "POS ");
		GetBlk(&m_file, &blk);
		SeekToType(&m_file,5, sizeof(float)*3);

		size=nbh*sizeof(float);
		P=new strParticleData [nbh];
		my_fread(P, size, 3, &m_file);
		swap_Nbyte((char*)P,nbh*3,4);
		/////////////////////////////////////////////
		cout<<"VEL"<<endl;
		sizeall=find_block(&m_file, "VEL ");
		GetBlk(&m_file, &blk);
		SeekToType(&m_file,4, sizeof(float)*3);	
		size=nbh*sizeof(float);
		pVel= new strParticleData[nbh];
		my_fread(pVel, size, 3, &m_file);
		swap_Nbyte((char*)pVel,nbh*3,4);
		/////////////////////////////////////////////
		int ibh=0;
		for(ip=0;ip<nbh;ip++)
			{
			m_data.push_back(
				strSPHParticle(&P[ip].Pos[0],char(5)));		
			cout<<"BHparticle: "<<ip<<") "<<P[ip].Pos[0]<<" "<<P[ip].Pos[1]<<" "<<P[ip].Pos[2]<<endl;
			m_data[ibh].sph.EnDt=0.0f;
			m_data[ibh].sph.Hsml=0.1f;
			m_data[ibh].sph.Rho=0.0f;
			m_data[ibh].sph.Temp=1.0f;

			ibh++;
			}
		m_NBH=nbh;
		/////////////////////////////////////////////
		

		
		}
	return m_isgood;
	}

///////////////////////////////////////////////////
// Read Stellar density if not exist then generate
bool CGadget::GetSPHParticles(int type, CRegion reg, bool flag_putin_COM)
	{
	unsigned int nstars=myhead.npart[type], ninreg=0;
	if(nstars<1) return false;
	unsigned int sizeall, size, ip;
	float *pM=NULL, *pH, *pF, *pU;
	strParticleData *pVel=NULL;
	cout<<"POS"<<endl;
	sizeall=find_block(&m_file, "POS ");
	GetBlk(&m_file, &blk);
	SeekToType(&m_file,type, sizeof(float)*3);

	size=nstars*sizeof(float);
	P=new strParticleData [nstars];
	my_fread(P, size, 3, &m_file);
	swap_Nbyte((char*)P,nstars*3,4);
	if(flag_putin_COM)
		PutInCOM(P, nstars);
	/////////////////////////////////////////////
	cout<<"VEL"<<endl;
	sizeall=find_block(&m_file, "VEL ");
	GetBlk(&m_file, &blk);
	SeekToType(&m_file,4, sizeof(float)*3);	
	size=nstars*sizeof(float);
	pVel= new strParticleData[nstars];
	my_fread(pVel, size, 3, &m_file);
	swap_Nbyte((char*)pVel,nstars*3,4);
	/////////////////////////////////////////////
	for(ip=0;ip<nstars;ip++)
		{
		if(reg.PointInRect(P[ip].Pos))
			{	
			m_data.push_back(
				strSPHParticle(&P[ip].Pos[0]));
			indexinreg.push_back(ip);			
			}
		}
	ninreg=m_data.size();	
	/////////////////////////////////////////////
	char buf[10];
	sprintf(buf,"%d", type);
	string rhofilename=m_filename+string("_rho_")+string(buf);
	bool get_rho_from_file=FileExists(rhofilename);
	ifstream rhofile;
	rhofile.open(rhofilename.c_str(),ios::in|ios::binary);
	bool o_swap=swp_flag;
	GetFileFormat(rhofile);
	//	if(get_rho_from_file)
	//{				
		  //io_header rhohead;
		  //GetHead(&rhofile, rhohead);
		//get_rho_from_file=(rhohead.npart[type]==ninreg);
	//cout<<swp_flag<<" before"<<o_swap<<endl;
	//exit(0);
	//		}
	
	if(get_rho_from_file){
/////////////////////////////////////////////
	  bool good_data=true;
		sizeall=find_block(&rhofile, rhoname[type].c_str());
		cout<< rhoname[type]<<endl;
		pF=new float[myhead.npart[type]];
		GetBlk(&rhofile, &blk);
		good_data=m_isgood;
	 	my_fread(pF, sizeof(float)*myhead.npart[type], 1, &rhofile);		
		swap_Nbyte((char*)pF,myhead.npart[type],4);

		pH=new float[myhead.npart[type]];
		sizeall=find_block(&rhofile,  hsmlname[type].c_str());
		cout<< hsmlname[type]<<endl;
		GetBlk(&rhofile, &blk);		
		good_data&=m_isgood;
		if(!good_data&&m_isgood){cout<<"Cannot get Rho block:"<<uname[type].c_str()<<endl;exit(0);};
	
		my_fread(pH, sizeof(float)*myhead.npart[type], 1, &rhofile);
		swap_Nbyte((char*)pH,myhead.npart[type],4);

		pU=new float[myhead.npart[type]];
		sizeall=find_block(&rhofile, uname[type].c_str());
		cout<<uname[type]<<endl;
		GetBlk(&rhofile, &blk);		
		my_fread(pU, sizeof(float)*myhead.npart[type], 1, &rhofile);
		swap_Nbyte((char*)pU,myhead.npart[type],4);
		if(!good_data&&m_isgood){cout<<"Cannot get Hsml block:"<<uname[type].c_str()<<
				"  for particle type "<<type<<endl;exit(0);};
		swp_flag=o_swap;
		for(unsigned i=0;i<ninreg;i++)
		  {
		    ip=indexinreg[i];
		    m_data[i].sph.Rho=	pF[ip];
		    m_data[i].sph.Hsml=	pH[ip];
		    m_data[i].sph.Temp=	pU[ip];
		  }
		delete [] pF;
		delete [] pH;
		delete [] pU;
		unit_conversion();
		/////////////////////////////////////////////
		}else{
#ifdef MYTREE			

			float one_particle_mass=1.0f;
			bool getmass_flag=false;
			if(myhead.mass[type]==0 && myhead.npart[type]>0)
				getmass_flag=true;
			else
				one_particle_mass=(float)myhead.mass[type];
			/////////////////////////////////////////////
			
			
			pM=new float[nstars];
			if(getmass_flag)
				{
				cout<<"MASS"<<endl;
				sizeall=find_block(&m_file, "MASS");
				GetBlk(&m_file, &blk);
				SeekToType(&m_file,4, sizeof(float), true);
				size=nstars*sizeof(float);
				my_fread(pM, size, 1, &m_file);		
				swap_Nbyte((char*)pM,nstars,4);
				}
			cout<<"Fill Sellar Matter"<<endl;
			for(unsigned int i=0;i<ninreg;i++)
				{
				if(getmass_flag)
					{	
					one_particle_mass=pM[indexinreg[i]];
					}	
				ip=indexinreg[i];
				v3 r(P[ip].Pos[0],P[ip].Pos[1],P[ip].Pos[2]), v(&pVel[ip].Pos[0]);			
				//v3 r(data[i].Pos[0],data[i].Pos[1],data[i].Pos[2]), v(&pVel[ip].Pos[0]);							
				float V2=pVel[ip].Pos[0]*pVel[ip].Pos[0];
				m_matter.push(CStone(r, v,one_particle_mass, i, V2));			
				}
			/////////////////////////////////////////////////
			//CheckRhoFile(reg);
			///////////////// Make Rho /////////////////////
			float vBox=reg.m_R*2;
			v3 center(reg.m_cent[0],reg.m_cent[1],reg.m_cent[2]);

			CBHTree *bhtree= new CBHTree(m_matter, vBox);

			//To test just generate some particles here   by GenFakeParticles();
			if(!m_matter.empty())
				{
				cout<<"Sending particles to tree: Ntotal= "<<m_matter.size()<<endl;
				CTree *pTree=new CTree(m_matter, center, vBox);	 
				cout<<"End tree"<<endl;
				cout<<"Making Rho by tree: Nnodes= "<<pTree->tree.size()<<endl;
				pTree->SetRSmooth(40);
				pTree->MakeRhoByTree();	 	 
				////////////////////////////////////////////////
				int count=0;
				for( list<node>::iterator no = pTree->tree.begin(); (no != pTree->tree.end()); ++no) 
					{
					if (no->leaf)
						{		
						//////////////////////////
						m_data[no->Type].sph.Rho=no->Rho;
						m_data[no->Type].sph.Temp=no->Temp;
						m_data[no->Type].sph.Hsml=no->radius;
						//cout<<no->Type<<")  "<<no->r<<" "<<no->Rho<<endl;
						count++;
						//////////////////////////			 
						}
					}
				cout<<"\nCount report sould be same:"<<count<<" "<<
					m_data.size()<<" Node Factor:"<<
					pTree->tree.size()/float(m_data.size())<<endl;
				delete pTree;
				}else{
					//make rho by BHTree;					
					for(unsigned int i=0;i<ninreg;i++)
						{
						m_data[i].sph.Rho =(float)bhtree->pp[(ninreg-1)-i].rho;
						m_data[i].sph.Temp=(float)bhtree->pp[(ninreg-1)-i].temp;
						m_data[i].sph.Hsml=(float)bhtree->pp[(ninreg-1)-i].hsml;
						}
					};
			if(pM!=NULL)delete [] pM;	
			//WriteRhoFile(rhofilename, type);
#endif

		}

	delete [] P;
	P=NULL;


	indexinreg.clear();
	return true;
	};
/*ReadOneBlock*/
unsigned int CGadget::read_block(float *&pV, char *name, int t)
{
  pV=new float[myhead.npart[t]];
  unsigned int sizeall=find_block(&m_file, name);
  cout<<name<<endl;
  if(sizeall<1) return 0;
  GetBlk(&m_file, &blk);		
  my_fread(pV, sizeof(float)*myhead.npart[t], 1, &m_file);
  swap_Nbyte((char*)pV,myhead.npart[t],4);
  return myhead.npart[t];

};
/*ReadOneBlock*/
unsigned int CGadget::read_blockv3(float *&pV, char *name, int t)
{
  pV=new float[myhead.npart[t]*3];
  unsigned int sizeall=find_block(&m_file, name);
  cout<<name<<endl;
  if(sizeall<1) return 0;
  GetBlk(&m_file, &blk);		
  SeekToType(&m_file,t, sizeof(float)*3);
  my_fread(pV, sizeof(float)*myhead.npart[t], 3, &m_file);
  swap_Nbyte((char*)pV,myhead.npart[t]*3,4);
  return myhead.npart[t];

};

/*Write One block*/
void CGadget::WriteOneBlock(ostream &file,string blname, char* pData, unsigned int datasize)
	{
		unsigned int blsize, idata;
		/*write block name*/
		blsize=8;
		file.write((char*)&blsize,4);
		file.write(blname.c_str(),4);		
		idata=8+datasize;
		file.write((char*)&idata,4);
		file.write((char*)&blsize,4);
/////////////////////////////////
		/*Writing data*/
		blsize=datasize;
		file.write((char*)&blsize,4);
		file.write((char*)pData,datasize);
		file.write((char*)&blsize,4);	 

	}
/* Write Stellar Rho File into _rho file*/
void CGadget::WriteRhoFile(string rhofilename, int type)
	{
		ofstream file;
		file.open(rhofilename.c_str(),ios::binary);
		unsigned int blsize, i,ninreg=m_data.size();
		io_header head;
		float *pHsml=new float[ninreg];
		float *pRho=new float[ninreg];
		float *pU=new float[ninreg];
		file.clear();
		memset(&head,0,sizeof(head));
		memset(head.npart,0,sizeof(head.npart));
		memset(head.npartTotal,0,sizeof(head.npart));
		head.npart[type]=ninreg;
		head.npartTotal[type]=ninreg;
		head.num_files=1;
/////////////////////////////////////////////////////
		blsize=sizeof(head);
		WriteOneBlock(file,"HEAD", (char*)&head, blsize);
//////////////////////////////////////////////////////
		for(i=0;i<ninreg;i++)
			{
			pRho[i]=m_data[i].sph.Rho;
			pHsml[i]=m_data[i].sph.Hsml;	
			pU[i]=m_data[i].sph.Temp;				
			}
//////////////////////////////////////////////////////		
		blsize=sizeof(float)*ninreg;
		WriteOneBlock(file,string("RHO "), (char*)pRho, blsize);
		WriteOneBlock(file,string("HSML"), (char*)pHsml, blsize);		
		WriteOneBlock(file,string("U   "), (char*)pU, blsize);
/////////////////////////////////////////////////////
		file.close();
		delete [] pRho;
		delete [] pHsml;
		delete [] pU;
	}
/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int  CGadget::unit_conversion(void)
{
  //    double GRAVITY, BOLTZMANN, PROTONMASS;
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
    double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
    double G, Xh, HubbleParam;

    unsigned int i;
    double MeanWeight, u, gamma;
    double RhoMean=1.9*1.0e-29;// in g*cm^-3;
    double OmegaBar=0.04;
  
	/* physical constants in cgs units */
    double    GRAVITY   = 6.672e-8,
        BOLTZMANN = 1.3806e-16,
        PROTONMASS = 1.6726e-24;

    /* internal unit system of the code */
    UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
    UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
    UnitVelocity_in_cm_per_s= 1.0e5;

    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
    UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
    UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

    G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);


    Xh= 0.76;  /* mass fraction of hydrogen */
    HubbleParam= 0.65;


    for(i=0; i<m_data.size(); i++)
    {
#ifdef SFR
	MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;
#else
	    MeanWeight= 4.0/(3*Xh+1+4*Xh) * PROTONMASS;
#endif

            /* convert internal energy to cgs units */

            u  = m_data[i].sph.Temp * 
				UnitEnergy_in_cgs/ UnitMass_in_g;

            gamma= 5.0/3;

            /* get temperature in Kelvin */

            m_data[i].sph.Temp= (float)
				(MeanWeight/BOLTZMANN * (gamma-1) * u);
			m_data[i].sph.Rho= float(
				m_data[i].sph.Rho* UnitDensity_in_cgs);

			m_data[i].sph.Rho/=float(RhoMean*OmegaBar);
    }
	return 0;
}


bool CGadget::ReadData(string file)
	{
	m_file.open(file.data(),  ios::in|ios::binary);

	char name[5];
	memset(name, 0,sizeof(name));
	if(m_file.bad())
		return false;

	
	GetBlk(&m_file, &blk);
	GetBlkName(&m_file, name);



	GetBlk(&m_file, &blk);
	GetHeader(&m_file);

	int nid=0, id=0, sizeall, size, nskip=0;
	GetBlk(&m_file, &nid);
	sizeall=find_block(&m_file, "ID  ");
	size=this->myhead.npart[4]*sizeof(int);

	if(sizeall !=size)
		{
		SeekToType(&m_file,4, sizeof(int));
		}
	ID=new int [this->myhead.npart[4]];
	GetBlk(&m_file, &blk);
	my_fread(ID, size, 1, &m_file);
	swap_Nbyte((char*)ID,this->myhead.npart[4],4);
	
	sizeall=find_block(&m_file, "POS ");
	GetBlk(&m_file, &blk);
	size=this->myhead.npart[4]*sizeof(float);
	P=new strParticleData [this->myhead.npart[4]];

	SeekToType(&m_file,4, sizeof(strParticleData));
	my_fread(P, size, 3, &m_file);
	swap_Nbyte((char*)P,this->myhead.npart[4]*3,4);


	return true;
	}
void CGadget::GetHeader(ifstream *fd){
	m_name=string("HEAD");
	find_block(fd,m_name.c_str());
	GetBlk(fd, &blk);
	my_fread((void*)myhead.npart,6*sizeof(int), 1, fd);             swap_Nbyte((char*)myhead.npart,6,4);
	my_fread((void*)myhead.mass,6*sizeof(double), 1, fd);           swap_Nbyte((char*)myhead.mass,6,8);
	my_fread((void*)&myhead.time,sizeof(double), 1, fd);            swap_Nbyte((char*)&myhead.time,1,8);
	my_fread((void*)&myhead.redshift,sizeof(double), 1, fd);        swap_Nbyte((char*)&myhead.redshift,1,8);
	my_fread((void*)&myhead.flag_sfr,sizeof(int), 1, fd);           swap_Nbyte((char*)&myhead.flag_sfr,1,4);
	my_fread((void*)&myhead.flag_feedback,sizeof(int), 1, fd);      swap_Nbyte((char*)&myhead.flag_feedback,1,4);
	my_fread((void*)myhead.npartTotal,6*sizeof(int), 1, fd);        swap_Nbyte((char*)myhead.npartTotal,6,4);
	my_fread((void*)&myhead.flag_cooling,sizeof(int), 1, fd);       swap_Nbyte((char*)&myhead.flag_cooling,1,4);
	my_fread((void*)&myhead.num_files,sizeof(int), 1, fd);          swap_Nbyte((char*)&myhead.num_files,1,4);
	my_fread((void*)&myhead.BoxSize,sizeof(double), 1, fd);         swap_Nbyte((char*)&myhead.BoxSize,1,8);
	my_fread((void*)&myhead.Omega0,sizeof(double), 1, fd);          swap_Nbyte((char*)&myhead.Omega0,1,8);
	my_fread((void*)&myhead.OmegaLambda,sizeof(double), 1, fd);     swap_Nbyte((char*)&myhead.OmegaLambda,1,8);
	my_fread((void*)&myhead.HubbleParam,sizeof(double), 1, fd);     swap_Nbyte((char*)&myhead.HubbleParam,1,8);
	my_fread((void*)&myhead.flag_multiphase,sizeof(int), 1, fd);    swap_Nbyte((char*)&myhead.flag_multiphase,1,4);
	my_fread((void*)&myhead.flag_stellarage,sizeof(int), 1, fd);    swap_Nbyte((char*)&myhead.flag_stellarage,1,4);
	my_fread((void*)&myhead.flag_sfrhistogram,sizeof(int), 1, fd);  swap_Nbyte((char*)&myhead.flag_sfrhistogram,1,4);
	my_fread((void*)&myhead.flag_metals,sizeof(int), 1, fd);        swap_Nbyte((char*)&myhead.flag_metals,1,4);
	my_fread((void*)&myhead.flag_decouple,sizeof(int), 1, fd);      swap_Nbyte((char*)&myhead.flag_decouple,1,4);
	my_fread((void*)&myhead.flag_effmodel,sizeof(int), 1, fd);      swap_Nbyte((char*)&myhead.flag_effmodel,1,4);
	my_fread((void*)myhead.fill,72*sizeof(char), 1, fd);
	GetBlk(fd, &blk);
	cout<<"====================="<<endl;
	for(unsigned it=0;it<6;it++)
	{
		printf("N[%d]=%0.9d\tMass[%d]=%g\n",it, myhead.npart[it], it, 
				myhead.mass[it]);
	};
	cout<<"====================="<<endl;
	};
void CGadget::GetHead(ifstream *fd, io_header &head){
	m_name=string("HEAD");
	find_block(fd,m_name.c_str());
	GetBlk(fd, &blk);
	my_fread((void*)head.npart,6*sizeof(int), 1, fd);             swap_Nbyte((char*)head.npart,6,4);
	my_fread((void*)head.mass,6*sizeof(double), 1, fd);           swap_Nbyte((char*)head.mass,6,8);
	my_fread((void*)&head.time,sizeof(double), 1, fd);            swap_Nbyte((char*)&head.time,1,8);
	my_fread((void*)&head.redshift,sizeof(double), 1, fd);        swap_Nbyte((char*)&head.redshift,1,8);
	my_fread((void*)&head.flag_sfr,sizeof(int), 1, fd);           swap_Nbyte((char*)&head.flag_sfr,1,4);
	my_fread((void*)&head.flag_feedback,sizeof(int), 1, fd);      swap_Nbyte((char*)&head.flag_feedback,1,4);
	my_fread((void*)head.npartTotal,6*sizeof(int), 1, fd);        swap_Nbyte((char*)head.npartTotal,6,4);
	my_fread((void*)&head.flag_cooling,sizeof(int), 1, fd);       swap_Nbyte((char*)&head.flag_cooling,1,4);
	my_fread((void*)&head.num_files,sizeof(int), 1, fd);          swap_Nbyte((char*)&head.num_files,1,4);
	my_fread((void*)&head.BoxSize,sizeof(double), 1, fd);         swap_Nbyte((char*)&head.BoxSize,1,8);
	my_fread((void*)&head.Omega0,sizeof(double), 1, fd);          swap_Nbyte((char*)&head.Omega0,1,8);
	my_fread((void*)&head.OmegaLambda,sizeof(double), 1, fd);     swap_Nbyte((char*)&head.OmegaLambda,1,8);
	my_fread((void*)&head.HubbleParam,sizeof(double), 1, fd);     swap_Nbyte((char*)&head.HubbleParam,1,8);
	my_fread((void*)&head.flag_multiphase,sizeof(int), 1, fd);    swap_Nbyte((char*)&head.flag_multiphase,1,4);
	my_fread((void*)&head.flag_stellarage,sizeof(int), 1, fd);    swap_Nbyte((char*)&head.flag_stellarage,1,4);
	my_fread((void*)&head.flag_sfrhistogram,sizeof(int), 1, fd);  swap_Nbyte((char*)&head.flag_sfrhistogram,1,4);
	my_fread((void*)&head.flag_metals,sizeof(int), 1, fd);        swap_Nbyte((char*)&head.flag_metals,1,4);
	my_fread((void*)&head.flag_decouple,sizeof(int), 1, fd);      swap_Nbyte((char*)&head.flag_decouple,1,4);
	my_fread((void*)&head.flag_effmodel,sizeof(int), 1, fd);      swap_Nbyte((char*)&head.flag_effmodel,1,4);
	my_fread((void*)head.fill,72*sizeof(char), 1, fd);
	GetBlk(fd, &blk);

	};

bool CGadget::GetFileFormat(ifstream &filein)
	{
	  swp_flag=false;
	GetBlk(&filein, &blk);
	if(blk != 8)
		{
		swp_flag=true;
		swap_Nbyte((char*)&blk,1,4);
		if(blk!=8)
			{
			cout<<"Cannot get file format..."<<endl;
			return false;
			}
		}
	filein.seekg(0, ios_base::beg);
	return true;	
	}
bool CGadget::GetFileFormat()
	{
	GetBlk(&m_file, &blk);
	if(blk != 8)
		{
		swp_flag=true;
		swap_Nbyte((char*)&blk,1,4);
		if(blk!=8)
			{
			cout<<"Cannot get file format..."<<endl;
			return false;
			}
		}
	m_file.seekg(0, ios_base::beg);
	return true;	
	}

//////////////////////////////////////////

//////////////////////////////////////
//////////////////////////////////////////
// Comparators for CGalaxy objects
// Ascending date sorting function
class cmpADMassSort
	{
	public:
		cmpADMassSort(bool mode=true):m_mode(mode){};
		bool operator()( CGalaxy*const& rpStart,  CGalaxy*const& rpEnd)
			{
			if(m_mode)
				return (rpStart->Ntotal) < (rpEnd->Ntotal);
			else
				return (rpStart->Ntotal) > (rpEnd->Ntotal);
			}
	private: 
		bool m_mode;
	};

// Descending data sorting function
struct cmpDescendingMassSort
	{
	bool operator()(CGalaxy*& rpStart, CGalaxy*& rpEnd)
		{
		return rpStart->Ntotal > rpEnd->Ntotal;
		}
	};
//////////////////////////////////////

struct HalogrpId: public std::binary_function< CGalaxy, int, bool > {
	bool operator () ( const CGalaxy &gal, const int &grpID ) const {
		return gal.grpID == grpID;
		}
	};


// this is our function object class.
// it holds a reference to a vector of grpIDS.
class CgrpIDScomp : public std::binary_function<unsigned int,unsigned int,bool> 
	{
	// the vector reference
	const vector<unsigned int>	&m_grpIDS;

	public:
		// constructor which takes a reference to a vector of CgrpIDS
		CgrpIDScomp( const vector<unsigned int> & CgrpIDS ) : m_grpIDS(CgrpIDS) {}

		// comparison operator. this will be called by std::sort with two
		// integers. we'll perform a comparison operation on the inputs and
		// return the result of the comparison.
		bool operator()(int a, int b) const
			{
			// a typical comparison operator compares the two input parameters.
			// this one uses the two input parameters as indexes into the m_grpIDS vector.
			//			return (m_grpIDS.at(a).ID) < (m_grpIDS.at(b).ID);

			return (m_grpIDS.at(a)) < (m_grpIDS.at(b));
			}
	};



CFOFCatalog::CFOFCatalog(string file):m_afofcat(file)
	{
	cout<<"Reading catalog file :"<<m_afofcat<<endl;
	if(!(m_isgood=ReadCatalog(m_afofcat)))
		{
		cout<<"Cannot open catalog file :"<<m_afofcat<<endl;
		cout<<"Exiting..."<<endl;
		exit(14);
		}


	}
bool CFOFCatalog::ReadCatalogIDsBin(string filename)
	{
	unsigned int nid=0, i=0;

	ifstream ifile(filename.c_str(),ios::in|ios::binary);
	if(	ifile.fail())
		{
		string mes=string("Can not open file: " );
		mes+= filename;
		cout<<mes.c_str()<<endl;
		return false;
		}

	ifile.read((char*)(&nid), sizeof(unsigned int));
	grpIDs.resize(nid);
	cout<<"..";
	ifile.read((char*)(&grpIDs[0]), sizeof(unsigned int)*nid);
	ifile.close();
	indexVector.clear();
	for(i=0;i<nid;i++)
		{
		// initialize indexVector
		if(grpIDs[i] >0)
			{
			indexVector.push_back(i);
			}
		}

	timer.start();
	cout<<".S.";
	sort(indexVector.begin(), 
		indexVector.end(), 
		CgrpIDScomp(grpIDs)); 
	timer.stop(DEBUG_ME);
	FillIDs();
	return true;
	};


class CPosDist2comp : public std::binary_function<int,int,bool> 
	{
	const vector<float>	&m_grpIDS;
	public:
		CPosDist2comp( const vector<float> & CgrpIDS ) : m_grpIDS(CgrpIDS) {}
		bool operator()(int a, int b) const
			{
			return (m_grpIDS.at(a)) < (m_grpIDS.at(b));
			}
	};

void CFOFCatalog::FillIDs()
	{
	unsigned int iindex=0, ind=0, ii=0;
	static bool first_flag=true;
	vector<CGalaxy*>::iterator it;// used for accessing to sorted elements
	cout<<".FID. "<<indexVector.size()<<"  ";
#ifdef DumpIDs
	FILE *testf;

	if(first_flag)
		{
		testf=fopen("C:/mingw/1.0/home/arm2arm/DATA/pos.txt", "w");
		first_flag=false;
		}else
			testf=fopen("C:/mingw/1.0/home/arm2arm/DATA/pos.txt", "a");

#endif
		vector<float> dist2;
		vector<int> indexOfDist2;
		float d2;
		for(iindex=0;iindex<indexVector.size();iindex++)
			{
			int value=grpIDs[indexVector[iindex]];

			for(it=halos.begin();it<halos.end();it++)
				{
				if((*it)->grpID == value)
					{
					(*it)->id.resize((*it)->Ntotal);
					dist2.clear();
					indexOfDist2.clear();
					for(ind=0;ind<(unsigned int)((*it)->Ntotal);ind++)
						{
						(*it)->id[ind]=IDs[indexVector[iindex+ind]];
						strParticleData *p=&(this->m_gadobj->P[indexVector[iindex+ind]]);
						indexOfDist2.push_back(ind);
						d2=0;
						for(int i=0;i<3;i++)				
							d2+= ((*p).Pos[i]-(*it)->x[i])*((*p).Pos[i]-(*it)->x[i]);//Sum(x-xc)^2
						dist2.push_back(d2);

						}
					sort(indexOfDist2.begin(), indexOfDist2.end(), CPosDist2comp(dist2));
					int iRpos=int(dist2.size()*RadFrac);
					(*it)->m_Rxx=sqrt(dist2.at(iRpos));
#ifdef DumpIDs
					if((*it)->ID==1)
						{
						fprintf(testf,"#%d  \n", iRpos);
						for(ind=0;ind<(unsigned int)((*it)->Ntotal);ind++)
							{
							strParticleData *p=&(this->m_gadobj->P[ indexVector[iindex+indexOfDist2[ind]] ] );
							fprintf(testf,"%f %f %f \n", (*p).Pos[0],(*p).Pos[1],(*p).Pos[2]);

							};
						}
#endif
					// Copy vector truncated by Rxx 
					(*it)->id_rxx.resize(iRpos);
					for(ind=0;ind<(unsigned int)(iRpos);ind++)
						{
						(*it)->id_rxx[ind]=(*it)->id[indexOfDist2[ind]];
						}
					//sort and check  
					(*it)->SortIds();
					iindex+=(*it)->Ntotal-1;
					//cout<<ii++<<"("<<(*it)->Ntotal<<")"<<"  ";
					}

				}
			}
#ifdef DumpIDs	
		fprintf(testf,"\n\n");
		fclose(testf);
#endif

	}
bool CFOFCatalog::ReadCatalogIDs(string file)
	{

	ifstream file_to_read(file.data());
	//	const int max_num_of_char_in_a_line = 512,// Maximum number of characters expected in a single line in the header
	//	num_of_header_lines = 1; // Number of header files to skip

	timer.start();
	if(DEBUG_ME)cout<<"Reading ID file: "<<file;
	if(file_to_read.bad())
		return false;
	int nid=0, id=0,i;

	file_to_read>>nid;
	///////////////////////////////
	grpIDs.resize(nid);
	///////////////////////////////
	for(i=0;i<nid;i++)
		{
		file_to_read>>grpIDs[i];
		// initialize indexVector
		if(grpIDs[i] >0)
			{
			indexVector.push_back(i);
			}
		}
	///////////////////////////////
	file_to_read.close();
	if(DEBUG_ME)
		cout<<"\ntotal "<<indexVector.size()<<"  and greather than zero: "<<float(indexVector.size())/float(nid)*100<<"perc... ok"<<endl;
	///////////////////////////////
	timer.stop(DEBUG_ME);
	// Do index sort by grpIDs
	if(DEBUG_ME)
		{
		int iindex=0;
		for(iindex=0;iindex<10;iindex++)
			cout<<indexVector[iindex]<<" "<<grpIDs[indexVector[iindex]]<<endl;


		cout<<"make heap"<<endl;
		timer.start();	
		//make_heap(indexVector.begin(), indexVector.end(),CgrpIDScomp(grpIDs));
		timer.stop(DEBUG_ME);		
		cout<<"sorting heap"<<endl;
		}
	timer.start();

	sort(indexVector.begin(), 
		indexVector.end(), 
		CgrpIDScomp(grpIDs)); 
	timer.stop(DEBUG_ME);
	FillIDs();
	return true;


	}

bool CFOFCatalog::ReadCatalog(string file)
	{
	ifstream file_to_read(file.data());
	const int max_num_of_char_in_a_line = 512,// Maximum number of characters expected in a single line in the header
		num_of_header_lines = 1; // Number of header files to skip
	if(DEBUG_ME)
		cout<<"Reading catalog file: "<<file;
	if(file_to_read.bad())
		return false;

	///////////////////////////////
	int Nhalo=0;
	file_to_read.ignore(1, '#');
	file_to_read>>Nhalo;
	file_to_read.ignore(max_num_of_char_in_a_line, '\n');
	CGalaxy* g;
	int tint=0;
	for( int ih=0;ih<Nhalo;ih++)
		{
		g=new CGalaxy;

		file_to_read>>(g->ID);
		file_to_read>>(g->grpID);
		file_to_read>>(g->x[0])>>(g->x[1])>>(g->x[2]);

		file_to_read>>(g->R);
		file_to_read>>(g->Ntotal);
		file_to_read>>(g->npart[0])>>(g->npart[1])>>(g->npart[2]);
		//		   file_to_read>>g->Rxyz[0]>>g->Rxyz[1]>>g->Rxyz[2];
		file_to_read.ignore(max_num_of_char_in_a_line, '\n');
		if(!file_to_read.good())return false;
		halos.push_back(g);
		}
	///////////////////////////////

	file_to_read.close();
	if(DEBUG_ME)
		cout<<"... ok"<<endl;

	if(DEBUG_ME)
		{
		cout<<"before sort"<<endl;
		for(int i=0;i<10;i++)
			cout<<(halos[i])->Ntotal<<endl;
		}
	if(DEBUG_ME)
		cout<<"Sorting catalog...";

	sort(halos.begin(), halos.end(),cmpADMassSort(false));

	if(DEBUG_ME)
		cout<<"...ok"<<endl;

	if(DEBUG_ME)
		{
		cout<<"after sort"<<endl;	
		for(int i=0;i<10;i++)
			cout<<(halos[i])->Ntotal<<endl;
		}

	return true;
	}















