#include "data_readers.h"

#define DEBUG_ME false
#define RadFrac  0.6
/////////////////////////////
CTimer timer; //timer for tests
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

bool StringToInt(const string &s, int &i) {
    istringstream myStream(s);

    if (myStream >> i)
        return true;
    else
        return false;
}

bool GetSnapPath(string snap, string &path) {
    basic_string <char>::size_type indexCh;
    static const basic_string <char>::size_type npos = 0;

    indexCh = snap.find_last_of('/');
    if ((indexCh) != npos)
        path.assign(snap, 0, int(indexCh));
    else
        path.assign("./");
    return true;
}

bool GetSnapName(string &snap, int &isnap) {
    string dig_snap;
    basic_string <char>::size_type indexCh;
    static const basic_string <char>::size_type npos = 0;

    indexCh = snap.find_last_of('/');
    if ((indexCh) != npos)
        snap.assign(snap, int(indexCh + 1), int(snap.size() - 1));

    while (snap.find_first_of("_") != 0) {
        indexCh = snap.find_first_of("_");
        dig_snap = dig_snap.assign(snap, int(indexCh + 1), int(4));
        if (StringToInt(dig_snap, isnap)) {
            snap = dig_snap;
            return true;
        } else
            snap.assign(snap, int(indexCh + 1), int(snap.size() - 1));

    }
    return false;
}
//////////////////////

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream) {
    size_t nread;

    if ((nread = fread(ptr, size, nmemb, stream)) != nmemb) {
        printf("I/O error (fread) %d !\n", (unsigned int) nread);
        exit(3);
    }
    return nread;
}

void CGadget::SeekToType(ifstream *file_to_read, int type, int onesize, bool mass_flag) {
    long nskip = 0;
    for (int ip = 0; ip < type; ip++) {
        if (!mass_flag)
            nskip += (this->myhead.npart[ip]);
        else {
            bool flag = (myhead.npart[ip] * myhead.mass[ip]) == 0.0 ? true : false;
            nskip += (this->myhead.npart[ip] * flag);
        }
    }

    file_to_read->seekg((nskip) * onesize, ios_base::cur);
    //cout<<"Iwill skip:"<<nskip<<endl;
    //exit(0);
}

int CGadget::find_block(ifstream *fd, const char *label) {
    int4bytes blocksize = 0, blksize = 0, ret = 0;
    char blocklabel[5] = {"    "};
    find_flag = true;
    fd->seekg(ios::beg);
    printf("Finding: %s\n", label);
    while (!fd->eof() && (blocksize == 0))//&& strcmp(blocklabel, "Z   ")!=0 )
    {
        //      cout<<"SKIP1"<<endl;
        GetBlk(fd, &blksize);
        //      cout<<"SKIP2"<<endl;
        //				if(blksize == 134217728)
        {
            swap_Nbyte((char*) &blksize, 1, 4);
        }
        if (blksize != 8) {
            if (ret > 0) {
                printf("incorrect format (blksize=%d)!\n", blksize);
                exit(1);
            } else break;
        } else {

            ret = my_fread(blocklabel, 4 * sizeof (char), 1, fd);
            if (ret > 0)
                ret = my_fread(&blocksize, sizeof (int4bytes), 1, fd);
            else
                return false;
            swap_Nbyte((char*) &blocksize, 1, 4);
            if (DEBUG_ME) {
                if (ret > 0)
                    printf("Found Block <%s> with %d bytes\n", blocklabel, blocksize);
                else
                    printf(" <%s> \n", label);


            }

            GetBlk(fd, &blksize);
            if (strcmp(label, blocklabel) != 0) {
                fd->seekg(blocksize, ios_base::cur);
                blocksize = 0;
            }
        }
    }
    find_flag = false;
    return (blocksize - 8);
}

bool CGadget::ReinitMultiFileName(int f) {
    m_file.close();
    sprintf(buffer, ".%d", f);
    nfilename = m_filename + string(buffer);
    cout << "In file: " << nfilename;
    m_file.open(nfilename.c_str(), ios::in | ios::binary);
    GetFileFormat();
    GetHeader();
    return true; // TODO make some sanity checks.
}
// Warning this function is not Thread Safe

bool CGadget::GetOneFile(CRegion reg, bool flag_putin_COM, bool isID) {
    unsigned int ngas = myhead.npart[0], ninreg = 0;
    unsigned int sizeall, size, ip;
    float *pF;
    CRange range, range1;
    sizeall = find_block(&m_file, "POS ");
    cout << "POS" << endl;
    GetBlk(&m_file, &blk);
    size = ngas * sizeof (float);
    P = new strParticleData [ngas];
    cout << "got" << ngas << endl;
    my_fread(P, size, 3, &m_file);
    swap_Nbyte((char*) P, ngas * 3, 4);
    if (flag_putin_COM) {
        PutInCOM(P, ngas);
    }
    range.Reset();
    for (ip = 0; ip < ngas; ip++) {
        range.getbound(P[ip].Pos[0]);
        range.getbound(P[ip].Pos[1]);
        range.getbound(P[ip].Pos[2]);
    }
    range.print((char*) "Pos Bounds: ");
    unsigned int data_shift = m_data.size();
    for (ip = 0; ip < ngas; ip++) {
        if (reg.PointInRect(P[ip].Pos)&& (m_frac > mydrand48())) {
            indexinreg.push_back(ip);
            m_data.push_back(
                    strSPHParticle(&P[ip].Pos[0]));
            ninreg++;
        }
    }

    // New2014 This is unused for the moment 
    // but you can read ID of particles for the post processing
    // This may help to read multiple files with non gas particles
    if (isID) {
        /////////////////////////////////////////////
        cout << "ID" << endl;
        sizeall = find_block(&m_file, "ID  ");
        GetBlk(&m_file, &blk);
        size = ngas * sizeof (int);
        int *pID = new int[ngas];
        my_fread(pID, size, 1, &m_file);
        swap_Nbyte((char*) pID, ngas, 4);
        /////////////////////////////////////////////
        for (ip = 0; ip < ninreg; ip++) {
            m_data[ip + data_shift].sph.ID = pID[indexinreg[ip]];
            //cout<<m_data[ip+data_shift].sph.ID<<endl;
        }
        delete [] pID;
    }
    // End New2014


    sizeall = find_block(&m_file, "RHO ");
    cout << "RHO" << endl;
    GetBlk(&m_file, &blk);
    my_fread(P, size, 1, &m_file);
    pF = (float *) P;
    swap_Nbyte((char*) pF, ngas, 4);
    for (ip = 0; ip < ninreg; ip++) {
        m_data[ip + data_shift].sph.Rho = pF[indexinreg[ip]];
        m_data[ip + data_shift].sph.Temp = 0.0f;
    }
    sizeall = find_block(&m_file, "U   ");
    cout << "U" << endl;
    GetBlk(&m_file, &blk);
    my_fread(P, size, 1, &m_file);
    pF = (float *) P;
    swap_Nbyte((char*) pF, ngas, 4);
    for (ip = 0; ip < ninreg; ip++) {
        m_data[ip + data_shift].sph.Temp = pF[indexinreg[ip]];
    }

    sizeall = find_block(&m_file, "HSML");
    cout << "HSML" << endl;
    GetBlk(&m_file, &blk);
    my_fread(P, size, 1, &m_file);
    pF = (float *) P;
    swap_Nbyte((char*) pF, ngas, 4);
    for (ip = 0; ip < ninreg; ip++) {
        m_data[ip + data_shift].sph.Hsml = pF[indexinreg[ip]];
    }
    /*
        sizeall = find_block(&m_file, "ENDT");
        cout << "ENDT" << endl;
        GetBlk(&m_file, &blk);
        my_fread(P, size, 1, &m_file);
        pF = (float *) P;
        swap_Nbyte((char*) pF, ngas, 4);
        for (ip = 0; ip < ninreg; ip++) {
            m_data[ip+data_shift].sph.EnDt = pF[indexinreg[ip]];
            range.getbound(m_data[ip+data_shift].sph.EnDt);
        }
        range.print((char*) "Entropy range:");
     */

    delete [] P;
    P = NULL;
    indexinreg.clear();
    return true;

}

bool CGadget::GetGas(CRegion reg, bool flag_putin_COM) {
    m_data.clear();
    if (m_multifile) {
        for (int f = 0; f < myhead.num_files; f++) {
            ReinitMultiFileName(f);
            GetOneFile(reg, flag_putin_COM);
        }

    } else {
        GetOneFile(reg, flag_putin_COM);
    }

    cout << "Total particles in region: " << m_data.size() << endl;

    CRange range, range1, range2;
    unsigned int ip = 0, ninreg = m_data.size();
    range.Reset();
    range1.Reset();
    range2.Reset();
    unit_conversion();
    for (ip = 0; ip < ninreg; ip++) {
        range.getbound(m_data[ip].sph.Rho);
        range1.getbound(m_data[ip].sph.Temp);
        range2.getbound(m_data[ip].sph.Hsml);
    }
    range.print((char*) "Rho  range:");
    range1.print((char*) "Temp range:");
    range2.print((char*) "Hsml range:");

};

// Check rhoi file if exist then read else make rho
// the rho file is inside the file snapfile+"_rho"
//WARNING!!!doesnot work with multiple file format.

bool CGadget::CheckRhoFile(CRegion reg) {
    string rhofilename = m_filename;
    if (m_multifile) {
        string rhofilename = m_filename + "_rho";
    }
    if (!FileExists(rhofilename)) {
        cout << "Can not find Rho file to get density: \n" << rhofilename << endl;

        return false;
    };

    return true;
}

bool CGadget::GetStars(CRegion reg, bool flag_putin_COM) {

    return GetSPHParticles(4, reg, flag_putin_COM);
}

bool CGadget::GetMYGas(CRegion reg, int Type) {

    return GetSPHParticles(Type, reg, false); //flag_putin_COM);	
}

bool CGadget::GetBH(CRegion reg) {
    if (this->myhead.npart[5] > 10)
        ; //GetSPHParticles(5, reg,false);	
    else {
        unsigned int nbh = myhead.npart[5];
        unsigned int sizeall, size, ip;
        float *pM = NULL;
        strParticleData *pVel = NULL;
        cout << "POS" << endl;
        sizeall = find_block(&m_file, "POS ");
        GetBlk(&m_file, &blk);
        SeekToType(&m_file, 5, sizeof (float)*3);

        size = nbh * sizeof (float);
        P = new strParticleData [nbh];
        my_fread(P, size, 3, &m_file);
        swap_Nbyte((char*) P, nbh * 3, 4);
        /////////////////////////////////////////////
        cout << "VEL" << endl;
        sizeall = find_block(&m_file, "VEL ");
        GetBlk(&m_file, &blk);
        SeekToType(&m_file, 4, sizeof (float)*3);
        size = nbh * sizeof (float);
        pVel = new strParticleData[nbh];
        my_fread(pVel, size, 3, &m_file);
        swap_Nbyte((char*) pVel, nbh * 3, 4);
        /////////////////////////////////////////////

        for (ip = 0; ip < nbh; ip++) {

            m_data.push_back(
                    strSPHParticle(&P[ip].Pos[0], char(5)));

            m_data.back().sph.EnDt = 0.0f;
            m_data.back().sph.Hsml = 0.2f;
            m_data.back().sph.Rho = 0.0f;
            m_data.back().sph.Temp = 0.0f;
            cout << "Got a BH particle: " << m_data.back() << endl;
        }
        cout << "Total BH particles: " << nbh << endl;
        m_NBH = nbh;
        /////////////////////////////////////////////
        delete [] pVel;


    }
    return m_isgood;
}


//Warning This is not thread safe function!!!

bool CGadget::ReinitMultiFileNameRho(int f, int Type) {
    char tbuffer[100];
    m_file.close();
    sprintf(tbuffer, "_rho_%d.%d", Type, f);
    nfilename = m_filename + string(tbuffer);
    cout << "In file: " << nfilename << endl;

    m_file.open(nfilename.c_str(), ios::in | ios::binary);

    GetFileFormat(m_file);

    GetHeader();

    return true; // TODO make some sanity checks.
}

///////////////////////////////////////////////////
// Read Stellar density if not exist then generate(TODO)
// TODO 2015 Add multiple file read function

bool CGadget::GetSPHParticles(int type, CRegion reg, bool flag_putin_COM) {

    if (type == 0) {
        cerr << "ERROR in data_readers.cpp:CGadget::GetSPHParticles:\n Please do not use this function for normal GAS!!!" << endl;
        exit(0);
    }

    unsigned int npart = 0, ninreg = 0;

    unsigned int sizeall, size, ip;

    int *TypeArr;
    //strParticleData *pVel = NULL;
    //The algorithm:
    //1) read pos in reg with IDps in reg
    //2) read Rho/Hsml from separate files selecting by IDps in reg


    // m_data.clear();
    cout << "POS" << endl;

    if (m_multifile) {
        for (int f = 0; f < myhead.num_files; f++) {
            ReinitMultiFileName(f);
            //ReinitMultiFileNameRho(f, type);
            //GetPosFromFile(reg, flag_putin_COM);


            sizeall = find_block(&m_file, "POS ");
            GetBlk(&m_file, &blk);
            if (type < 6) {
                npart = myhead.npart[type];
            } else {
                for (int it = 0; it < 6; it++)
                    npart += myhead.npart[it];
            }

            if (npart < 1) continue;

            TypeArr = new int[npart];

            if (type < 6)
                std::fill(TypeArr, TypeArr + npart, type);
            else {
                int ip = 0;
                for (int it = 0; it < 6; it++) {
                    for (int i = 0; i < myhead.npart[it]; i++) {
                        TypeArr[ip] = it;
                        ip++;
                    }
                }
            }


            if (type < 6)
                SeekToType(&m_file, type, sizeof (float)*3);

            size = npart * sizeof (float);
            P = new strParticleData [npart];

            my_fread(P, size, 3, &m_file);
            swap_Nbyte((char*) P, npart * 3, 4);


            /////////////////////////////////////////////
            cout << "ID" << endl;
            sizeall = find_block(&m_file, "ID  ");
            GetBlk(&m_file, &blk);
            if (type < 6)
                SeekToType(&m_file, type, sizeof (int));
            size = npart * sizeof (int);
            int *pID = new int[npart];
            my_fread(pID, size, 1, &m_file);

            //////////////////////////////////////////////
            //Get Velocity
            /*
            sizeall = find_block(&m_file, "VEL ");
            GetBlk(&m_file, &blk);
            if (type < 6)
                SeekToType(&m_file, type, sizeof (float)*3);
            size = npart * sizeof (float);
            strParticleData *pVel = new strParticleData [npart];
            my_fread(pVel, size, 3, &m_file);
            swap_Nbyte((char*) pVel, npart * 3, 4);
            */
            //////////////////////////////////////////////
            cout <<"In CGadget::GetSPHParticles: "<< reg.m_cent[0] << " ";
            cout << reg.m_cent[1] << " ";
            cout << reg.m_cent[2] << endl;
            //          double vcom[]={ -432.62475, 291.90301,-197.61080};
            double vcom[] = {0.0, 0.0, 0.0};
            /////////////////////////////////////////////            
            for (ip = 0; ip < npart; ip++) {

                if (reg.PointInRect(P[ip].Pos)) {
                    m_data.push_back(
                            strSPHParticle(&P[ip].Pos[0]));
                    m_data[m_data.size() - 1].sph.ID = pID[ip];
                   /* pVel[ip].Pos[0] -= vcom[0];
                    pVel[ip].Pos[1] -= vcom[1];
                    pVel[ip].Pos[2] -= vcom[2];
                    m_data[m_data.size() - 1].sph.Temp = sqrt(pVel[ip].Pos[0] * pVel[ip].Pos[0] +
                            pVel[ip].Pos[1] * pVel[ip].Pos[1] +
                            pVel[ip].Pos[2] * pVel[ip].Pos[2]);
                    m_data[m_data.size() - 1].sph.Rho = m_data[m_data.size() - 1].sph.Temp;
*/
                    m_data[m_data.size() - 1].sph.Type = type;
                    idinreg[ pID[ip] ] = m_data.size() - 1;
                }

            }//end ip loop 
            delete [] pID;
            delete [] P;
            //delete [] pVel;
            P = NULL;
            delete [] TypeArr;
        }//end of files loop 

    }// end is multifile format if

    ninreg = m_data.size();
    cout << "npart " << ninreg << endl;


    /////////////////////////////////////////////
    if (ninreg > 1) {
        //if (false) {    

        long ichanged = 0;
        for (int f = 0; f < myhead.num_files; f++) {
            //	ReinitMultiFileName(f);
            bool isgood = ReinitMultiFileNameRho(f, type);
            //TODO:
            //Probably we need to check NP in file and NP rhofile 
            npart = myhead.npart[0];
            cout << "Type=" << type << " " << rhoname[type] << " " << npart << endl;

            int *pID = new int[npart];
            float *pH = NULL, *pF = NULL;
            pF = new float[npart];
            pH = new float[npart];

            sizeall = find_block(&m_file, "ID  ");
            GetBlk(&m_file, &blk);
            my_fread(pID, sizeall, 1, &m_file);
            //cout << "Sall=" << sizeall << " inreg=" << idinreg.size() << endl;

            sizeall = find_block(&m_file, rhoname[type].c_str());
            GetBlk(&m_file, &blk);
            my_fread(pF, sizeall, 1, &m_file);

            sizeall = find_block(&m_file, hsmlname[type].c_str());
            GetBlk(&m_file, &blk);
            my_fread(pH, sizeall, 1, &m_file);

            for (unsigned i = 0; i < npart; i++) {
                //cout <<"rh="<<pF[i] <<"  h="<<pH[i]<< endl;
                if (idinreg.count(pID[i]) > 0) {
                    ip = idinreg.at(pID[i]);
                    m_data[ip].sph.Hsml = pH[i];
                    m_data[ip].sph.Rho = pF[i];
                    //m_data[ip].sph.Temp = pF[i];

                    //m_data[i].sph.Type = TypeArr[ip]; 
                    // cout <<"rh="<<pF[i] <<"  h="<<pH[i]<< endl;
                    ichanged++;
                }

            }

            delete [] pF;
            delete [] pH;
            delete [] pID;

            //unit_conversion();
        }
        if (ninreg != ichanged) {
            cout << "Something is odd with RHO file:" << endl;
            cout << "npart " << ninreg << " but changed" << ichanged << endl;
            exit(1);
        }
    }
#ifdef PPS    
    /////////////////////////////////////////////
    //Some Post processing of data.....
    printStats();
#define GRID 200

    ////////////////////////////////////////////////////////
    typedef vector<unsigned int> iarr;
    std::vector<vector<vector<iarr> > > igrid;

    // Set up sizes. (HEIGHT x WIDTH)
    igrid.resize(GRID);
    for (int i = 0; i < GRID; ++i) {
        igrid[i].resize(GRID);
        for (int j = 0; j < GRID; ++j)
            igrid[i][j].resize(GRID);
    }
printStats();

    ////////////////////////////////////////////////////////    
    unsigned int ix[3], ii = 0;
    for (int i = 0; i < m_data.size(); i++) {
        for (ii = 0; ii < 3; ii++) {
            ix[ii] = (m_data[i].Pos[ii] - vrange[ii].Min) / (vrange[ii].Max - vrange[ii].Min) * GRID;            
            if (ix[ii] < 0)ix[ii] = 0;
            if (ix[ii] > GRID-1)ix[ii] = GRID - 1;
            
        }
        igrid[ix[0]][ix[1]][ix[2]].push_back(i);
    }

    for (int i = 0; i < m_data.size(); i++) {
        float h = m_data[i].sph.Hsml;
        int cx = 0, cy, cz;
        for (int ix = -h; ix < h; ix++)
            for (int iy = -h; iy < h; iy++)
                for (int iz = -h; iz < h; iz++) {
                    if (cx < ix && cy < iy && cz < iz) {
                        m_data[i].sph.Temp;
                    }

                }
    }

    exit(0);
#endif    
    return true;
};

/*ReadOneBlock*/
unsigned int CGadget::read_block(float *&pV, const char *name, int t)
{
 
  
  int nall=myhead.npart[t];
  
  if(t==6)
    {
        nall=0;
        for(int i=0;i<6;i++)
            nall+=myhead.npart[i];
    }else 
      if(t=!0)nall=myhead.npart[0];
  pV=new float[nall];
  unsigned int sizeall=find_block(&m_file, name);
  cout<<name<<endl;
  if(sizeall<1) return 0;
  GetBlk(&m_file, &blk);		
  my_fread(pV, sizeof(float)*nall, 1, &m_file);
  swap_Nbyte((char*)pV,nall,4);
  return nall;

};
/*ReadOneBlock*/
unsigned int CGadget::read_blockv3(float *&pV, const char *name, int t)
{
    int nall=myhead.npart[t];    
    if(t==6)
    {
        nall=0;
        for(int i=0;i<6;i++)
            nall+=myhead.npart[i];
    }
  pV=new float[nall*3];
  unsigned int sizeall=find_block(&m_file, name);
  cout<<name<<endl;
  if(sizeall<1) return 0;
  GetBlk(&m_file, &blk);
  if(t !=6)SeekToType(&m_file,t, sizeof(float)*3);
  else SeekToType(&m_file,0, sizeof(float)*3);
  my_fread(pV, sizeof(float)*nall, 3, &m_file);
  swap_Nbyte((char*)pV,nall*3,4);
  return nall;

};
/*Write One block*/
void CGadget::WriteOneBlock(ostream &file, string blname, char* pData, unsigned int datasize) {

    unsigned int blsize, idata;
    /*write block name*/
    blsize = 8;
    file.write((char*) &blsize, 4);
    file.write(blname.c_str(), 4);
    idata = 8 + datasize;
    file.write((char*) &idata, 4);
    file.write((char*) &blsize, 4);
    /////////////////////////////////
    /*Writing data*/
    blsize = datasize;
    file.write((char*) &blsize, 4);
    file.write((char*) pData, datasize);
    file.write((char*) &blsize, 4);

}

/* Write Stellar Rho File into _rho file*/
void CGadget::WriteRhoFile(string rhofilename, int type) {
    ofstream file;
    file.open(rhofilename.c_str(), ios::binary);
    unsigned int blsize, i, ninreg = m_data.size();
    io_header head;
    float *pHsml = new float[ninreg];
    float *pRho = new float[ninreg];
    float *pU = new float[ninreg];
    file.clear();
    memset(&head, 0, sizeof (head));
    memset(head.npart, 0, sizeof (head.npart));
    memset(head.npartTotal, 0, sizeof (head.npart));
    head.npart[type] = ninreg;
    head.npartTotal[type] = ninreg;
    head.num_files = 1;
    /////////////////////////////////////////////////////
    blsize = sizeof (head);
    WriteOneBlock(file, "HEAD", (char*) &head, blsize);
    //////////////////////////////////////////////////////
    for (i = 0; i < ninreg; i++) {

        pRho[i] = m_data[i].sph.Rho;
        pHsml[i] = m_data[i].sph.Hsml;
        pU[i] = m_data[i].sph.Temp;
    }
    //////////////////////////////////////////////////////		
    blsize = sizeof (float)*ninreg;
    WriteOneBlock(file, string("RHO "), (char*) pRho, blsize);
    WriteOneBlock(file, string("HSML"), (char*) pHsml, blsize);
    WriteOneBlock(file, string("U   "), (char*) pU, blsize);
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
int CGadget::unit_conversion(void) {
    //    double GRAVITY, BOLTZMANN, PROTONMASS;
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
    double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
    double G, Xh, HubbleParam;

    unsigned int i;
    double MeanWeight, u, gamma;
    double RhoMean = 1.9 * 1.0e-29; // in g*cm^-3;
    double OmegaBar = 0.04;

    /* physical constants in cgs units */
    double GRAVITY = 6.672e-8,
            BOLTZMANN = 1.3806e-16,
            PROTONMASS = 1.6726e-24;

    /* internal unit system of the code */
    UnitLength_in_cm = 3.085678e21; /*  code length unit in cm/h */
    UnitMass_in_g = 1.989e43; /*  code mass unit in g/h */
    UnitVelocity_in_cm_per_s = 1.0e5;

    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
    UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
    UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

    G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);


    Xh = 0.76; /* mass fraction of hydrogen */
    HubbleParam = 0.65;


    for (i = 0; i < m_data.size(); i++) {
#ifdef SFR
        MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh * P[i].Ne) * PROTONMASS;
#else

        MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh) * PROTONMASS;
#endif

        /* convert internal energy to cgs units */

        u = m_data[i].sph.Temp *
                UnitEnergy_in_cgs / UnitMass_in_g;

        gamma = 5.0 / 3;

        /* get temperature in Kelvin */

        m_data[i].sph.Temp = (float)
                (MeanWeight / BOLTZMANN * (gamma - 1) * u);
        m_data[i].sph.Rho = float(
                m_data[i].sph.Rho * UnitDensity_in_cgs);

        m_data[i].sph.Rho /= float(RhoMean * OmegaBar);
    }
    return 0;
}

bool CGadget::ReadData(string file) {
    m_file.open(file.data(), ios::in | ios::binary);

    char name[5];
    memset(name, 0, sizeof (name));
    if (m_file.bad())
        return false;


    GetBlk(&m_file, &blk);
    GetBlkName(&m_file, name);



    GetBlk(&m_file, &blk);
    GetHeader(&m_file);

    int nid = 0, id = 0, sizeall, size, nskip = 0;
    GetBlk(&m_file, &nid);
    sizeall = find_block(&m_file, "ID  ");
    size = this->myhead.npart[4] * sizeof (int);

    if (sizeall != size) {
        SeekToType(&m_file, 4, sizeof (int));
    }
    ID = new int [this->myhead.npart[4]];
    GetBlk(&m_file, &blk);
    my_fread(ID, size, 1, &m_file);
    swap_Nbyte((char*) ID, this->myhead.npart[4], 4);

    sizeall = find_block(&m_file, "POS ");
    GetBlk(&m_file, &blk);
    size = this->myhead.npart[4] * sizeof (float);
    P = new strParticleData [this->myhead.npart[4]];

    SeekToType(&m_file, 4, sizeof (strParticleData));
    my_fread(P, size, 3, &m_file);
    swap_Nbyte((char*) P, this->myhead.npart[4]*3, 4);

    return true;
}

void CGadget::GetHeader(ifstream *fd) {
    m_name = string("HEAD");
    find_block(fd, m_name.c_str());
    GetBlk(fd, &blk);
    my_fread((void*) myhead.npart, 6 * sizeof (int), 1, fd);
    swap_Nbyte((char*) myhead.npart, 6, 4);
    my_fread((void*) myhead.mass, 6 * sizeof (double), 1, fd);
    swap_Nbyte((char*) myhead.mass, 6, 8);
    my_fread((void*) &myhead.time, sizeof (double), 1, fd);
    swap_Nbyte((char*) &myhead.time, 1, 8);
    my_fread((void*) &myhead.redshift, sizeof (double), 1, fd);
    swap_Nbyte((char*) &myhead.redshift, 1, 8);
    my_fread((void*) &myhead.flag_sfr, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_sfr, 1, 4);
    my_fread((void*) &myhead.flag_feedback, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_feedback, 1, 4);
    my_fread((void*) myhead.npartTotal, 6 * sizeof (int), 1, fd);
    swap_Nbyte((char*) myhead.npartTotal, 6, 4);
    my_fread((void*) &myhead.flag_cooling, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_cooling, 1, 4);
    my_fread((void*) &myhead.num_files, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.num_files, 1, 4);
    my_fread((void*) &myhead.BoxSize, sizeof (double), 1, fd);
    swap_Nbyte((char*) &myhead.BoxSize, 1, 8);
    my_fread((void*) &myhead.Omega0, sizeof (double), 1, fd);
    swap_Nbyte((char*) &myhead.Omega0, 1, 8);
    my_fread((void*) &myhead.OmegaLambda, sizeof (double), 1, fd);
    swap_Nbyte((char*) &myhead.OmegaLambda, 1, 8);
    my_fread((void*) &myhead.HubbleParam, sizeof (double), 1, fd);
    swap_Nbyte((char*) &myhead.HubbleParam, 1, 8);
    my_fread((void*) &myhead.flag_multiphase, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_multiphase, 1, 4);
    my_fread((void*) &myhead.flag_stellarage, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_stellarage, 1, 4);
    my_fread((void*) &myhead.flag_sfrhistogram, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_sfrhistogram, 1, 4);
    my_fread((void*) &myhead.flag_metals, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_metals, 1, 4);
    my_fread((void*) &myhead.flag_decouple, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_decouple, 1, 4);
    my_fread((void*) &myhead.flag_effmodel, sizeof (int), 1, fd);
    swap_Nbyte((char*) &myhead.flag_effmodel, 1, 4);
    my_fread((void*) myhead.fill, 72 * sizeof (char), 1, fd);
    GetBlk(fd, &blk);
    cout << "=====================" << endl;
    for (unsigned it = 0; it < 6; it++) {

        printf("N[%d]=%0.9d\tMass[%d]=%g\n", it, myhead.npart[it], it,
                myhead.mass[it]);
    };
    cout << "=====================" << endl;
};

void CGadget::GetHead(ifstream *fd, io_header &head) {

    m_name = string("HEAD");
    find_block(fd, m_name.c_str());
    GetBlk(fd, &blk);
    my_fread((void*) head.npart, 6 * sizeof (int), 1, fd);
    swap_Nbyte((char*) head.npart, 6, 4);
    my_fread((void*) head.mass, 6 * sizeof (double), 1, fd);
    swap_Nbyte((char*) head.mass, 6, 8);
    my_fread((void*) &head.time, sizeof (double), 1, fd);
    swap_Nbyte((char*) &head.time, 1, 8);
    my_fread((void*) &head.redshift, sizeof (double), 1, fd);
    swap_Nbyte((char*) &head.redshift, 1, 8);
    my_fread((void*) &head.flag_sfr, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_sfr, 1, 4);
    my_fread((void*) &head.flag_feedback, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_feedback, 1, 4);
    my_fread((void*) head.npartTotal, 6 * sizeof (int), 1, fd);
    swap_Nbyte((char*) head.npartTotal, 6, 4);
    my_fread((void*) &head.flag_cooling, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_cooling, 1, 4);
    my_fread((void*) &head.num_files, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.num_files, 1, 4);
    my_fread((void*) &head.BoxSize, sizeof (double), 1, fd);
    swap_Nbyte((char*) &head.BoxSize, 1, 8);
    my_fread((void*) &head.Omega0, sizeof (double), 1, fd);
    swap_Nbyte((char*) &head.Omega0, 1, 8);
    my_fread((void*) &head.OmegaLambda, sizeof (double), 1, fd);
    swap_Nbyte((char*) &head.OmegaLambda, 1, 8);
    my_fread((void*) &head.HubbleParam, sizeof (double), 1, fd);
    swap_Nbyte((char*) &head.HubbleParam, 1, 8);
    my_fread((void*) &head.flag_multiphase, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_multiphase, 1, 4);
    my_fread((void*) &head.flag_stellarage, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_stellarage, 1, 4);
    my_fread((void*) &head.flag_sfrhistogram, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_sfrhistogram, 1, 4);
    my_fread((void*) &head.flag_metals, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_metals, 1, 4);
    my_fread((void*) &head.flag_decouple, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_decouple, 1, 4);
    my_fread((void*) &head.flag_effmodel, sizeof (int), 1, fd);
    swap_Nbyte((char*) &head.flag_effmodel, 1, 4);
    my_fread((void*) head.fill, 72 * sizeof (char), 1, fd);
    GetBlk(fd, &blk);

};

bool CGadget::GetFileFormat(ifstream &filein) {
    swp_flag = false;
    GetBlk(&filein, &blk);
    if (blk != 8) {
        swp_flag = true;
        swap_Nbyte((char*) &blk, 1, 4);
        if (blk != 8) {
            cout << "Cannot get file format..." << endl;
            return false;
        }
    }
    filein.seekg(0, ios_base::beg);

    return true;
}

bool CGadget::GetFileFormat() {
    GetBlk(&m_file, &blk);
    if (blk != 8) {
        swp_flag = true;
        swap_Nbyte((char*) &blk, 1, 4);
        if (blk != 8) {
            cout << "Cannot get file format..." << endl;
            return false;
        }
    }
    m_file.seekg(0, ios_base::beg);

    return true;
}

void CGadget::printStats() {

    if (!good()) {
        cerr << "no gadget data in container!!!" << endl;
        return;
    }
    for (int i = 0; i < 6; i++)vrange[i].Reset();

    for (unsigned int i = 0; i < m_data.size(); i++) {
        vrange[0].getbound(m_data[i].Pos[0]);
        vrange[1].getbound(m_data[i].Pos[1]);
        vrange[2].getbound(m_data[i].Pos[2]);
        vrange[3].getbound(m_data[i].sph.Rho);
        vrange[4].getbound(m_data[i].sph.Hsml);
        vrange[5].getbound(m_data[i].sph.Temp);
    }
    cout << "====== DEBUG=====" << endl;
    vrange[0].print((char*) "PosX range is: ");
    vrange[1].print((char*) "PosY range is: ");
    vrange[2].print((char*) "PosZ range is: ");
    vrange[3].print((char*) "Rho range is: ");
    vrange[4].print((char*) "Hsml range is: ");
    vrange[5].print((char*) "Temp range is: ");
    cout << "======+++++=====" << endl;
}
//////////////////////////////////////////

ostream& operator<<(ostream& output, const strSPHParticle& p) {
    output << "Type[" << p.sph.Type << "]-->"
            << "POS( "
            << p.Pos[0] << " , "
            << p.Pos[1] << " , "
            << p.Pos[2] << " ) <=> ";
    output << " hsml= " << p.sph.Hsml
            << " rho= " << p.sph.Rho
            << " temp= " << p.sph.Temp
            << " EnDt= " << p.sph.EnDt;

    return output; // for multiple << operators.
}


//////////////////////////////////////
//////////////////////////////////////////
// Comparators for CGalaxy objects
// Ascending date sorting function

class cmpADMassSort {
public:

    cmpADMassSort(bool mode = true) : m_mode(mode) {
    };

    bool operator()(CGalaxy * const& rpStart, CGalaxy * const& rpEnd) {
        if (m_mode)
            return (rpStart->Ntotal) < (rpEnd->Ntotal);

        else
            return (rpStart->Ntotal) > (rpEnd->Ntotal);
    }
private:
    bool m_mode;
};

// Descending data sorting function

struct cmpDescendingMassSort {

    bool operator()(CGalaxy*& rpStart, CGalaxy*& rpEnd) {

        return rpStart->Ntotal > rpEnd->Ntotal;
    }
};
//////////////////////////////////////

struct HalogrpId : public std::binary_function< CGalaxy, int, bool > {

    bool operator()(const CGalaxy &gal, const int &grpID) const {

        return gal.grpID == grpID;
    }
};


// this is our function object class.
// it holds a reference to a vector of grpIDS.

class CgrpIDScomp : public std::binary_function<unsigned int, unsigned int, bool> {
    // the vector reference

    const vector<unsigned int> &m_grpIDS;

public:
    // constructor which takes a reference to a vector of CgrpIDS

    CgrpIDScomp(const vector<unsigned int> & CgrpIDS) : m_grpIDS(CgrpIDS) {
    }

    // comparison operator. this will be called by std::sort with two
    // integers. we'll perform a comparison operation on the inputs and
    // return the result of the comparison.

    bool operator()(int a, int b) const {
        // a typical comparison operator compares the two input parameters.
        // this one uses the two input parameters as indexes into the m_grpIDS vector.
        //			return (m_grpIDS.at(a).ID) < (m_grpIDS.at(b).ID);

        return (m_grpIDS.at(a)) < (m_grpIDS.at(b));
    }
};

CFOFCatalog::CFOFCatalog(string file) : m_afofcat(file) {
    cout << "Reading catalog file :" << m_afofcat << endl;
    if (!(m_isgood = ReadCatalog(m_afofcat))) {

        cout << "Cannot open catalog file :" << m_afofcat << endl;
        cout << "Exiting..." << endl;
        exit(14);
    }


}

bool CFOFCatalog::ReadCatalogIDsBin(string filename) {
    unsigned int nid = 0, i = 0;

    ifstream ifile(filename.c_str(), ios::in | ios::binary);
    if (ifile.fail()) {
        string mes = string("Can not open file: ");
        mes += filename;
        cout << mes.c_str() << endl;
        return false;
    }

    ifile.read((char*) (&nid), sizeof (unsigned int));
    grpIDs.resize(nid);
    cout << "..";
    ifile.read((char*) (&grpIDs[0]), sizeof (unsigned int)*nid);
    ifile.close();
    indexVector.clear();
    for (i = 0; i < nid; i++) {
        // initialize indexVector
        if (grpIDs[i] > 0) {
            indexVector.push_back(i);
        }
    }

    timer.start();
    cout << ".S.";
    sort(indexVector.begin(),
            indexVector.end(),
            CgrpIDScomp(grpIDs));
    timer.stop(DEBUG_ME);
    FillIDs();

    return true;
};

class CPosDist2comp : public std::binary_function<int, int, bool> {
    const vector<float> &m_grpIDS;
public:

    CPosDist2comp(const vector<float> & CgrpIDS) : m_grpIDS(CgrpIDS) {
    }

    bool operator()(int a, int b) const {

        return (m_grpIDS.at(a)) < (m_grpIDS.at(b));
    }
};

void CFOFCatalog::FillIDs() {
    unsigned int iindex = 0, ind = 0, ii = 0;
    static bool first_flag = true;
    vector<CGalaxy*>::iterator it; // used for accessing to sorted elements
    cout << ".FID. " << indexVector.size() << "  ";
#ifdef DumpIDs
    FILE *testf;

    if (first_flag) {
        testf = fopen("C:/mingw/1.0/home/arm2arm/DATA/pos.txt", "w");
        first_flag = false;
    } else
        testf = fopen("C:/mingw/1.0/home/arm2arm/DATA/pos.txt", "a");

#endif
    vector<float> dist2;
    vector<int> indexOfDist2;
    float d2;
    for (iindex = 0; iindex < indexVector.size(); iindex++) {
        int value = grpIDs[indexVector[iindex]];

        for (it = halos.begin(); it < halos.end(); it++) {
            if ((*it)->grpID == value) {
                (*it)->id.resize((*it)->Ntotal);
                dist2.clear();
                indexOfDist2.clear();
                for (ind = 0; ind < (unsigned int) ((*it)->Ntotal); ind++) {
                    (*it)->id[ind] = IDs[indexVector[iindex + ind]];
                    strParticleData *p = &(this->m_gadobj->P[indexVector[iindex + ind]]);
                    indexOfDist2.push_back(ind);
                    d2 = 0;
                    for (int i = 0; i < 3; i++)
                        d2 += ((*p).Pos[i]-(*it)->x[i])*((*p).Pos[i]-(*it)->x[i]); //Sum(x-xc)^2
                    dist2.push_back(d2);

                }
                sort(indexOfDist2.begin(), indexOfDist2.end(), CPosDist2comp(dist2));
                int iRpos = int(dist2.size() * RadFrac);
                (*it)->m_Rxx = sqrt(dist2.at(iRpos));
#ifdef DumpIDs
                if ((*it)->ID == 1) {
                    fprintf(testf, "#%d  \n", iRpos);
                    for (ind = 0; ind < (unsigned int) ((*it)->Ntotal); ind++) {
                        strParticleData *p = &(this->m_gadobj->P[ indexVector[iindex + indexOfDist2[ind]] ]);
                        fprintf(testf, "%f %f %f \n", (*p).Pos[0], (*p).Pos[1], (*p).Pos[2]);

                    };
                }
#endif
                // Copy vector truncated by Rxx 
                (*it)->id_rxx.resize(iRpos);
                for (ind = 0; ind < (unsigned int) (iRpos); ind++) {
                    (
                            *it)->id_rxx[ind] = (*it)->id[indexOfDist2[ind]];
                }
                //sort and check  
                (*it)->SortIds();
                iindex += (*it)->Ntotal - 1;
                //cout<<ii++<<"("<<(*it)->Ntotal<<")"<<"  ";
            }

        }
    }
#ifdef DumpIDs	
    fprintf(testf, "\n\n");
    fclose(testf);
#endif

}

bool CFOFCatalog::ReadCatalogIDs(string file) {

    ifstream file_to_read(file.data());
    //	const int max_num_of_char_in_a_line = 512,// Maximum number of characters expected in a single line in the header
    //	num_of_header_lines = 1; // Number of header files to skip

    timer.start();
    if (DEBUG_ME)cout << "Reading ID file: " << file;
    if (file_to_read.bad())
        return false;
    int nid = 0, id = 0, i;

    file_to_read >> nid;
    ///////////////////////////////
    grpIDs.resize(nid);
    ///////////////////////////////
    for (i = 0; i < nid; i++) {
        file_to_read >> grpIDs[i];
        // initialize indexVector
        if (grpIDs[i] > 0) {
            indexVector.push_back(i);
        }
    }
    ///////////////////////////////
    file_to_read.close();
    if (DEBUG_ME)
        cout << "\ntotal " << indexVector.size() << "  and greather than zero: " << float(indexVector.size()) / float(nid)*100 << "perc... ok" << endl;
    ///////////////////////////////
    timer.stop(DEBUG_ME);
    // Do index sort by grpIDs
    if (DEBUG_ME) {
        int iindex = 0;
        for (iindex = 0; iindex < 10; iindex++)
            cout << indexVector[iindex] << " " << grpIDs[indexVector[iindex]] << endl;


        cout << "make heap" << endl;
        timer.start();
        //make_heap(indexVector.begin(), indexVector.end(),CgrpIDScomp(grpIDs));
        timer.stop(DEBUG_ME);
        cout << "sorting heap" << endl;
    }
    timer.start();

    sort(indexVector.begin(),
            indexVector.end(),
            CgrpIDScomp(grpIDs));
    timer.stop(DEBUG_ME);
    FillIDs();

    return true;


}

bool CFOFCatalog::ReadCatalog(string file) {
    ifstream file_to_read(file.data());
    const int max_num_of_char_in_a_line = 512, // Maximum number of characters expected in a single line in the header
            num_of_header_lines = 1; // Number of header files to skip
    if (DEBUG_ME)
        cout << "Reading catalog file: " << file;
    if (file_to_read.bad())
        return false;

    ///////////////////////////////
    int Nhalo = 0;
    file_to_read.ignore(1, '#');
    file_to_read >> Nhalo;
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    CGalaxy * g;
    int tint = 0;
    for (int ih = 0; ih < Nhalo; ih++) {
        g = new CGalaxy;

        file_to_read >> (g->ID);
        file_to_read >> (g->grpID);
        file_to_read >> (g->x[0]) >> (g->x[1]) >> (g->x[2]);

        file_to_read >> (g->R);
        file_to_read >> (g->Ntotal);
        file_to_read >> (g->npart[0]) >> (g->npart[1]) >> (g->npart[2]);
        //		   file_to_read>>g->Rxyz[0]>>g->Rxyz[1]>>g->Rxyz[2];
        file_to_read.ignore(max_num_of_char_in_a_line, '\n');
        if (!file_to_read.good())return false;
        halos.push_back(g);
    }
    ///////////////////////////////

    file_to_read.close();
    if (DEBUG_ME)
        cout << "... ok" << endl;

    if (DEBUG_ME) {
        cout << "before sort" << endl;
        for (int i = 0; i < 10; i++)
            cout << (halos[i])->Ntotal << endl;
    }
    if (DEBUG_ME)
        cout << "Sorting catalog...";

    sort(halos.begin(), halos.end(), cmpADMassSort(false));

    if (DEBUG_ME)
        cout << "...ok" << endl;

    if (DEBUG_ME) {
        cout << "after sort" << endl;
        for (int i = 0; i < 10; i++)
            cout << (halos[i])->Ntotal << endl;
    }

    return true;
}














