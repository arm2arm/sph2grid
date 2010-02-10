#ifdef  _OPENMP
#include <omp.h>
#endif

#include "sphvolume.h"
#include "data_readers.h"
#include "utils.h"
#include "Render.h"

bool allocated_flag = false;
float ***vol3d = NULL;
float ***vol3dsm = NULL;
int np, nx, ny, nz;
float *X, *Y, *Z, *hsml, *rho;

#define EXIT_FAIL exit(1);

void saveFree(float* v) {
    if (v != NULL) {
        delete [] v;
        v = NULL;
    }
}
void printStat(float *pV, const char* txt, int np) {
    CRange range;
    range.Reset();
    for (unsigned int ip = 0; ip < np; ip++) {
        range.getbound(pV[ip]);
    }
    range.print(txt);
}

float Wsph(float rr, float h) {
    float retval = 0;
    float u = rr / h;
    if (0 <= u && u <= 1)
        retval = 1.0f - 3.0f / 2.0f * u * u + 3.0f / 4.0f * u * u * u;
    else
        if (1 < u && u <= 2)
        retval = 1.0f / 4.0f * (2.0f - u)*(2.0f - u)*(2.0f - u);
    else
        return 0;

    return retval; //(3.1456f*h*h*h);


}

void DeallocateVol3D(float ***vol3d) {
    if (vol3d != NULL) {
        for (int i = 0; i < nz; ++i) {
            for (int j = 0; j < ny; ++j)
                delete [] vol3d[i][j];
            delete [] vol3d[i];
        }

        delete [] vol3d;
        vol3d = NULL;
    }
}

void ZeroVol3D(float ***&vol) {
    int i, j, k;
#ifdef  _OPENMP
#pragma omp parallel default(shared) 
#endif
    {
#ifdef  _OPENMP
#pragma omp for private(i,j,k)
#endif
        for (k = 0; k < nz; k++)
            for (i = 0; i < nx; i++)
                for (j = 0; j < ny; j++)
                    vol[i][j][k] = 0.0f;
    }
}

void CopyVol3D(float ***&src, float ***&dest) {
    int i, j, k;
#ifdef  _OPENMP
#pragma omp parallel default(shared) 
#endif
    {
#ifdef  _OPENMP
#pragma omp for private(i,j,k)
#endif
        for (k = 0; k < nz; k++)
            for (i = 0; i < nx; i++)
                for (j = 0; j < ny; j++)
                    dest[i][j][k] = src[i][j][k];
    }
}

void AllocateVol3D(float ***&vol) {
    int i, j, k;
    cout << "Allocate" << endl;
    if (vol != NULL) {
        cout << "Not null lets deallocate" << endl;
        DeallocateVol3D(vol);
    }
    if (!(nx && ny && nz)) {
        fprintf(stderr, "# Something is wrong with GRID size:"
                "\nGRIDX=%d\nGRIDY=%d\nGRIDZ=%d\n\n", nx, ny, nz);
        exit(13);
    }
    vol = new float** [nx];
    for (k = 0; k < nx; ++k) {
        vol[k] = new float* [ny];
        for (j = 0; j < ny; ++j)
            vol[k][j] = new float [nz];
    }

    cout << "End allocate" << endl;

}

int WriteSPHVolume(string fname) {
    ofstream fout;
    fout.open(fname.c_str(), ios::binary);
    if (!fout.fail()) {
        int ti[] = {nx, ny, nz};
        fout.write((char*) & ti[0], sizeof (int) *3);
        for (int k = 0; k < nz; k++) {
            float *slab = new float[nx * ny + 1];
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++) {
                    float trho = vol3d[i][j][k];
                    slab[i + j * nx] = trho;
                }

            fout.write((char*) & slab[0], sizeof (float) * nx * ny);
            // cout<<">>> Page: "<<k<<endl;
        }
        fout.close();

    } else
        cout << "  Cannot open the file: " << fname << endl;
    return 1;
}

int WriteSPHGRIDVolume(string fname, int type,
        float XC, float YC, float ZC, float RC, int GRID) {
    ofstream fout;
    fout.open(fname.c_str(), ios::binary);
    if (!fout.fail()) {
        int ti = GRID;
        float tf[] = {XC, YC, ZC, RC};
        fout.write((char*) & type, sizeof (int));
        fout.write((char*) & tf[0], sizeof (float) *4);
        fout.write((char*) & ti, sizeof (int));
        for (int k = 0; k < nz; k++) {
            float *slab = new float[nx * ny + 1];
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++) {
                    float trho = vol3d[i][j][k];
                    slab[i + j * nx] = trho;
                }

            fout.write((char*) & slab[0], sizeof (float) * nx * ny);
            // cout<<">>> Page: "<<k<<endl;
        }
        fout.close();
        return 0;
    } else
        cout << "  Cannot open the file: " << fname << endl;
    return 1;
}
///////////////////////

void DoSPHVolume() {
    float i, j, k;
    int ii, jj, kk;
    float r, h;
    ZeroVol3D(vol3d);
    cout << "starting SPHVol" << endl;
    for (int ip = 0; ip < np; ip++) {
        //      dist2=X[ip]*X[ip]+Y[ip]*Y[ip]+Z[ip]*Z[ip];
        if (ip % 10000 == 0)cout << ip / float(np)*100 << "%" << endl;
        h = hsml[ip];
        i = (X[ip]);
        j = (Y[ip]);
        k = (Z[ip]);
        if (i + h < nx && i - h > 0 && j + h < ny && j - h > 0 && k + h < nz && k - h > 0) {
#ifdef  _OPENMP
#pragma omp parallel default(shared) 
#endif
            {
#ifdef  _OPENMP
#pragma omp for private(kk, jj, ii,r)	 
#endif
                for (kk = int(k - h) + 1; kk<int(k + h) + 1; kk += 1) {
                    for (jj = int(j - h) + 1; jj<int(j + h) + 1; jj += 1)
                        for (ii = int(i - h) + 1; ii<int(i + h) + 1; ii += 1) {
                            r = sqrt((i - ii)*(i - ii)+(j - jj)*(j - jj)+(k - kk)*(k - kk));
                            vol3d[int(ii)][int(jj)][int(kk)] += rho[ip] * Wsph(r, h);
                        }

                }
                //if(ip%1000==0)cout<<ip<<" "<<h<<endl;
            }
        }
    }
}

void SmoothSPHVolume(int sm) {
    int i, j, k, ii, jj, kk;
    float rhomean;
    float count;
    int sm2 = int(0.5f * sm);
    int ks, js, is, ke, je, ie;

    ZeroVol3D(vol3dsm);
    ke = nz - sm2 - 1;
    je = ny - sm2 - 1;
    ie = nx - sm2 - 1;

    ks = sm2;
    js = sm2;
    is = sm2;

    cout << "Entering parallel region" << endl;
#ifdef  _OPENMP
#pragma omp parallel default(shared) 
#endif
    {
#ifdef  _OPENMP
#pragma omp master
#endif
        {
#ifdef _OPENMP 
            cout << "numthreads=" << omp_get_num_threads() << endl;
#endif
        }
#ifdef  _OPENMP
#pragma omp for private(kk, jj, ii,rhomean, count )
#endif
        for (kk = ks; kk < ke; kk++) {
            for (jj = js; jj < je; jj++)
                for (ii = is; ii < ie; ii++) {
                    rhomean = 0.0;
                    count = 0.0;
                    for (k = kk - sm2; k <= kk + sm2; k++)
                        for (j = jj - sm2; j <= jj + sm2; j++)
                            for (i = ii - sm2; i <= ii + sm2; i++) {
                                rhomean += vol3d[i][j][k];
                                count += 1.0;
                            }
                    //if(count>0)
                    vol3dsm[ii][jj][kk] = (rhomean / count);
                }

        };
    }
#ifdef  _OPENMP
#pragma omp parallel default(shared) 
#endif
    {
#ifdef  _OPENMP
#pragma omp for private(i,j,k)
#endif
        for (k = 0; k < nz; k++)
            for (i = 0; i < nx; i++)
                for (j = 0; j < ny; j++)
                    vol3d[i][j][k] += (vol3dsm[i][j][k]*0.5);
    }



}

void OldSmoothSPHVolume(int sm) {
    int i, j, k, ii, jj, kk;
    float rhomean;
    float count;
    int sm2 = int(0.5f * sm);
    int ks, js, is, ke, je, ie;

    ZeroVol3D(vol3dsm);
    CopyVol3D(vol3d, vol3dsm);
    ke = nz - sm2 - 1;
    je = ny - sm2 - 1;
    ie = nx - sm2 - 1;

    ks = sm2;
    js = sm2;
    is = sm2;

    for (kk = ks; kk < ke; kk++) {
        for (jj = js; jj < je; jj++)
            for (ii = is; ii < ie; ii++) {
                rhomean = 0.0;
                count = 0.0;
                for (k = kk - sm2; k <= kk + sm2; k++)
                    for (j = jj - sm2; j <= jj + sm2; j++)
                        for (i = ii - sm2; i <= ii + sm2; i++) {
                            rhomean += vol3d[i][j][k];
                            count += 1.0;
                        }

                vol3dsm[ii][jj][kk] += (rhomean / count);
            }

    }


    CopyVol3D(vol3dsm, vol3d);


}

void printusage(void) {
    cout << "***************************\n" <<
            "Usage: \n" <<
            "sph2vol.x Filename smoothingwindow Niter\n" <<
            "***************************\n\n" << endl;
}
///////////////

int rmain(int argc, char* argv[]) {

    ifstream filein;
    int smwindow = 2, NSmooth = 1;
    if (argc != 4) {
        cout << "\nError: incorrect number of parameters." << endl;
        printusage();
        exit(1);
    }

    filein.open(argv[1], ios::binary);


    filein.read((char*) (&np), sizeof (int));
    filein.read((char*) (&nx), sizeof (int));
    filein.read((char*) (&ny), sizeof (int));
    filein.read((char*) (&nz), sizeof (int));


    cout << "Np: " << np << " GRID= " << nx << " " << ny << " " << nz << endl;
    if (np > 1) {
        X = new float[np];
        Y = new float[np];
        Z = new float[np];
        rho = new float[np];
        hsml = new float[np];
        AllocateVol3D(vol3d);
        AllocateVol3D(vol3dsm);
    } else
        exit(13);
    /////////////////////////////////////////////////////
    filein.read((char*) (&X[0]), sizeof (float) * np);
    filein.read((char*) (&Y[0]), sizeof (float) * np);
    filein.read((char*) (&Z[0]), sizeof (float) * np);
    filein.read((char*) (&rho[0]), sizeof (float) * np);
    filein.read((char*) (&hsml[0]), sizeof (float) * np);
    filein.close();
    /////////////////////////////////////////////////////
    string fname = string(argv[1]);
    fname += string(".vol");
    DoSPHVolume();
    WriteSPHVolume(fname);

    smwindow = int(atof(argv[2]));
    cout << "Setting smoothing Window: " << smwindow << endl;
    if (smwindow <= 1)smwindow = 2;

    NSmooth = atoi(argv[3]);
    if (NSmooth <= 0)NSmooth = 1;


    for (int is = 0; is < NSmooth; is++) {
        cout << "starting smoothing SPHVol: " << is << endl;
#ifdef _OPENMP
        SmoothSPHVolume(smwindow);
#else
        OldSmoothSPHVolume(smwindow);
#endif
    }

    fname = string(argv[1]);
    fname += string(".sm.vol");
    cout << "writing...";
    WriteSPHVolume(fname);
    cout << "...done" << endl;
    delete [] X;
    delete [] Y;
    delete [] Z;
    delete [] hsml;
    delete [] rho;
    DeallocateVol3D(vol3d);
    DeallocateVol3D(vol3dsm);
    return 0;

}
/////////////////////////////

class CGetData {
public:

    CGetData(string file, int type = 0) : m_infile(file),
    m_type(type), pHSML(NULL),
    pRHO(NULL), pPOS(NULL),rho_loaded(false), hsml_loaded(false) {

        if (!ReadData()) {
            cout << "Fail to read data...exiting" << endl;
            //EXIT_FAIL;
        }
    };

    ~CGetData() {
        saveFree(pHSML);
        saveFree(pRHO);
        saveFree(pPOS);

    };

    unsigned int GetNp() {
        return nelem;
    };

    void SetAutoCOM(float &x, float &y, float &z) {
        float *pPot = NULL;
        CGadget *g = new CGadget(m_infile, false);
        int nparticles = g->read_block(pPot,
                "POT ",
                6);
        if (nparticles > 0) {
            int ip = 0;

            for (int i = 1; i < nparticles; i++)
                if (pPot[ip] > pPot[i]) {
                    ip = i;
                    //cout << i << " " << pPot[i] << endl;
                }
            saveFree(pPot);
            float *pP;

	    int nread=g->read_blockv3(pP,"POS ",6);

            if (nparticles != g->read_blockv3(pP,
                    "POS ",6)){
	      cout<<nread<<"!="<<nparticles<<endl;
	      EXIT_FAIL;
	    };
            std::cout <<"AutoCOM result:" <<nparticles << " " << ip
                    << " " << pP[ip*3]
                    << " " << pP[ip*3 + 1]
                    << " " << pP[ip*3 + 2] << endl;
           
            x = pP[ip*3];
            y = pP[ip*3 + 1];
            z = pP[ip*3 + 2];
            saveFree(pP);
           
            
        }

    };
public:
  bool rho_loaded, hsml_loaded;
protected:
 
    bool ReadData() {
        CGadget *m_gin = new CGadget(m_infile, false);
        bool ans;

        ans = GetPosBlock(m_gin);
        int npcurrent = GetNp(); //store the current number of particles
	if(m_type!=0)
	  {
	    string ext="_rho_"+intToString(m_type);
	    delete m_gin;
	    m_gin = new CGadget(m_infile+ext, false);
	    cout<<m_gin->m_isgood<<endl;
	  }
	if(m_gin->m_isgood)
	  {
	
	    if (!GetHsmlBlock(m_gin))nelem = npcurrent;
	    if (!GetRhoBlock(m_gin))nelem = npcurrent;
	    if(m_type==0)
	      if (!GetUBlock(m_gin))
		{
		  
		  nelem = npcurrent;
		}
	  }

        return ans;
    };

    bool GetRhoBlock(CGadget *g) {
        bool retval = true;

        nelem = 0;
        nelem = g->read_block(pRHO,
                (char *) g->rhoname[m_type].c_str(),
                m_type);

        if (nelem == 0) {
            cout << "Error reading block:" << g->rhoname[m_type].c_str() <<
                    " In file: " << g->m_filename << endl;
            retval = false;
        } else cout << "Getting: " << g->rhoname[m_type] << " with " << nelem << " elems.." << endl;
	rho_loaded=retval;
        return retval;

    };
   bool GetUBlock(CGadget *g) {
        bool retval = true;

        nelem = 0;
        nelem = g->read_block(pTEMP,
                "U   ",
                0);

        if (nelem == 0) {
            cout << "Error reading block: U " <<
                    " In file: " << g->m_filename << endl;
            retval = false;
        } else cout << "Getting: " << "U " << " with " << nelem << " elems.." << endl;

	double gamma= 5.0/3.0;
	double Xh= 0.76;  /* mass fraction of hydrogen */
	double BOLTZMANN = 1.3806e-16;
        double PROTONMASS = 1.6726e-24;
	double UnitMass_in_g= 1.989e43; 
	double UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
	double  UnitVelocity_in_cm_per_s= 1.0e5;
	double	UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  
	double   UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

	double MeanWeight= 4.0/(3*Xh+1+4*Xh) * PROTONMASS;
	float fac=(float)(MeanWeight/BOLTZMANN * (gamma-1.0)*
			  UnitEnergy_in_cgs/ UnitMass_in_g);
	for(int i=0;i<nelem;i++)
	  pTEMP[i]=(fac* pTEMP[i]);
        return retval;

    };

    bool GetHsmlBlock(CGadget *g) {
        bool retval = true;

        nelem = 0;
        nelem = g->read_block(pHSML,
                (char *) g->hsmlname[m_type].c_str(),
                m_type);
	//cout<<nelem<<endl;
	
        if (nelem == 0) {
            cout << "Error reading block:" << g->hsmlname[m_type].c_str() <<
                    " In file: " << g->m_filename << endl;
            retval = false;
        } else cout << "Getting: " << g->hsmlname[m_type] << " with " << nelem << " elems.." << endl;
	hsml_loaded=retval;
        return retval;

    };

    bool GetPosBlock(CGadget *g) {
        bool retval = true;

        nelem = 0;
        nelem = g->read_blockv3(pPOS,
                "POS ",
                m_type);

        if (nelem == 0) {
            cout << "Error reading block:" << "POS" <<
                    " In file: " << g->m_filename << endl;
            retval = false;
        } else cout << "Getting: " << "POS" << " with " << nelem << " elems.." << endl;
	printStat(pPOS, "pos range: ", nelem*3);
        return retval;

    };

protected:
    int m_type;
    unsigned int nelem;
    string m_infile;
public:
  float *pHSML;
  float *pTEMP;
  float *pRHO;
  float *pPOS;

};

void FreeSphMemory() {
    delete [] X;
    delete [] Y;
    delete [] Z;
    delete [] hsml;
    delete [] rho;
    DeallocateVol3D(vol3d);
    DeallocateVol3D(vol3dsm);
}

bool rhoSortFunction (int i,int j) { return (rho[i]>rho[j]); }

bool DoSph2Grid(string fname, string foutname, int type,
        float XC, float YC, float ZC, float RC, int GRID, float zfac = 1.0f, float hsml_in_kpc = 0.5) {
    bool volume_flag = false;
    bool render_flag = true;
    
    int *idx;

    /*Some checks */
    if (GRID < 16) {
        cout << "\nError:\nGrid should be > 16 you have: " << GRID << endl;
        EXIT_FAIL;
    }
    if (RC <= 0) {
        cout << "\nError:\nRc  should be > 0 you have:  " << RC << endl;
        EXIT_FAIL;
    }
    if (type != 0 && volume_flag) {
        cout << "\nError:\nIn this version of the code\n" <<
                "the Type  should be  0 you have:  " << type << endl;
        EXIT_FAIL;
    }
    /***********************/
    CGetData *data = new CGetData(fname, type);

    if (XC == 0 && YC == 0 && ZC == 0)
        data->SetAutoCOM(XC, YC, ZC);

    /***********************/
    np = data->GetNp();
    nx = ny = nz = GRID;
    nz =(int) (nz*zfac);

    cout << "Np: " << np << " GRID= " << nx << " " << ny << " " << nz << endl;
    if (np > 1) {
        X = new float[np];
        Y = new float[np];
        Z = new float[np];

        rho = new float[np];
        hsml = new float[np];
	idx = new int[np];
        if(volume_flag)
        AllocateVol3D(vol3d);
	std::cout<<"memory allocation done"<<std::endl;
    } else
        exit(13);

    if(data->rho_loaded)
      printStat(data->pRHO, "Rho range: ", np);
    if(data->hsml_loaded)
      printStat(data->pHSML, "HSML range: ", np);
    unsigned int ni = 0;
    float RC2 = RC * 2, BOX=64000.0f;
    float fixed_hsml = hsml_in_kpc * nx / RC; // Setup the smoothing length in kpc/h
    std::fill(X, X+np,0.0f);
    std::fill(Y, Y+np,0.0f);
    std::fill(Z, Z+np,0.0f);
    std::fill(rho, rho+np,1.0f);
    std::fill(hsml, hsml+np,std::max(1.0f,fixed_hsml));

    for (unsigned int i = 0, ic = 0; ic < np; ic++) {
        X[ni] = (data->pPOS[i ] - XC);
        Y[ni] = (data->pPOS[i + 1] - YC);
        Z[ni] = (data->pPOS[i + 2] - ZC);
       
        if (X[ni] < RC && X[ni] >-RC &&
                Y[ni] < RC && Y[ni] >-RC &&
                Z[ni] < RC*zfac && Z[ni] >-RC*zfac
                ) {
            X[ni] += RC;
            Y[ni] += RC;
            Z[ni] += RC;
            X[ni] *= nx / RC2;
            Y[ni] *= ny / RC2;
            Z[ni] *= nz / RC2;
	    //cout<<"  > "<<X[ni]<<endl;
	    if(data->rho_loaded && data->hsml_loaded)
	      {
		rho[ni] = data->pRHO[ic];///data->pTEMP[ic]*1e7;//*((nx / RC2)*(nx / RC2)*(nx / RC2));
		//rho[ni] = data->pTEMP[ic];//*((nx / RC2)*(nx / RC2)*(nx / RC2));
		hsml[ni] = std::max(1.0f,data->pHSML[ic] * nx / RC2);
	      }		
	    idx[ni]=ni;
            ni++;
        }
        i += 3;

    }
    
    std::sort(idx, idx+ni, rhoSortFunction);
    printStat(data->pPOS, "Pos range: ", ni);
    printStat(X, "X range: ", ni);
    printStat(Y, "Y range: ", ni);
    printStat(Z, "Z range: ", ni);
    printStat(hsml, "hsml range: ", ni);
    printStat(rho, "rho range: ", ni);
    cout << "Particles in region Ni=" << ni << " Out of total: " << np << endl;
   
    np = ni;
    if (volume_flag) {
        DoSPHVolume();

        cout << "Writing GRID data to the output file: " << foutname << endl;
        WriteSPHGRIDVolume(foutname, type,
                XC, YC, ZC, RC, GRID);
    }
    ////////////////////////////////////////////

    if (render_flag) {
        CRender *pRenderer = new CRender();
	if((int)getEnv("SPH2GRID_DUMPLUT",0))
	  {
	    pRenderer->DumpAllLUT();exit(0);
	  };
        pRenderer->SetColorTable((int)getEnv("SPH2GRID_LUT",14));
        if (volume_flag)
	  pRenderer->DoRenderByGrid(vol3d, GRID, zfac);
        else
	  {
	    pRenderer->SetMinMaxRho(
				    getEnv("SPH2GRID_MINRHO",-7.0f),
				    getEnv("SPH2GRID_MAXRHO",3.0f));
	    //	    pRenderer->DoRenderByPoints(X, Y, Z, std::max(fixed_hsml,1.0f), np, GRID);
	    std::cout<<"start to render"<<std::endl; 
	    if(getEnv("SPH2GRID_XYZ",1)){
	      pRenderer->SetOut("sph2grid.png");
	      if((int)getEnv("SPH2GRID_MAXINTENSITY",1))
		{
		  pRenderer->DoRenderByAdaptiveSortedPoints(X, Y,idx,rho, hsml, np, GRID);
		}else
		  pRenderer->DoRenderByAdaptivePoints(X, Y, Z, rho, hsml, np, GRID);
	      pRenderer->SetOut("sph2grid_90.png");
	      if((int)getEnv("SPH2GRID_MAXINTENSITY",1))
		pRenderer->DoRenderByAdaptiveSortedPoints(X, Z,idx,rho, hsml, np, GRID);
	      else
		pRenderer->DoRenderByAdaptivePoints(X, Z, Y, rho, hsml, np, GRID);
	    }else
	      {
		pRenderer->SetOut("sph2grid.png");
		if((int)getEnv("SPH2GRID_MAXINTENSITY",1))
		  pRenderer->DoRenderByAdaptiveSortedPoints(X, Y,idx,rho, hsml, np, GRID);
		else
		  pRenderer->DoRenderByAdaptivePoints(X, Y,Z, rho, hsml, np, GRID);
	      }
	  }
        delete pRenderer;
    }
    /////////////////////////////////////////////





    FreeSphMemory();
    delete data;


    return true;
}
