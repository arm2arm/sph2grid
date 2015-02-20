#ifndef _CRange_
#define _CRange_
#include <iostream>
#include <algorithm>
#include <string>
#include <cassert>
#include <vector>
using std::min;
using std::max;
using std::cout;
using std::endl;

class CRange{
public:
	float Min;
	float Max;
	float m_sum;
	float m_mean;
	unsigned int m_np;
	CRange(){Reset();};
	~CRange(){};
	void sum(float x){m_sum+=x;m_np++;};

	void getbound(float x)
	{
		Min=min(x, Min);
		Max=max(x, Max);
		sum(x);
	};
void Reset(void ){Min=1e10; Max=-1e10;m_sum=0;m_np=0;};
void print(const char *str){
	cout<<str<<" Min="<<Min<<"\tMax="<<Max<<"\t mean="<<m_sum/m_np<<endl;

};
};

std::string intToString(int i);
#define printOpenGLError() printOglError(__FILE__, __LINE__)
float getEnv(char * name, float defvalue);
const char * getEnvStr(char * name, const char defvalue[]);
int printOglError(char *file, int line);
float   stoptimer_( int *flag);
void  starttimer_( void);
float mydrand48(void);
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795f
#endif








template<class T>
class CKernel{
public:
        CKernel(const int dim=3):m_dim(dim){
                if((m_dim<1)||(m_dim>20))
                        {
                        std::cout<<"Specify the Normalization constant for dimensions > 20"<<endl;
                        
                        }
                AllocateTable();
                };
        ~CKernel(){};
        virtual T W(T u)=0;
protected:
        void AllocateTable()
                {
                KernelRad.resize(KERNEL_TABLE+2);//+2 allow to go in the loops  KERNEL_TABLE+1;
                Kernel.resize(KERNEL_TABLE+2);
                };
        virtual void    InitTable()=0;
protected:
        static const int KERNEL_TABLE=1000;
        std::vector<T> KernelRad,Kernel;
        int m_dim;
        std::string TypeOfKernel;
        };



template<class T>
class CEpanechikov : public CKernel<T>
        {
        public:
                CEpanechikov(int dim):CKernel<T>(dim){this->TypeOfKernel=std::string("Epanechikov kernel= Cd*(1-u)");InitTable();};
        private:
        virtual void    InitTable(){
                static double
                        fep[]={0.75000113,0.63661975,0.59683102,0.60792705,0.66492015,0.77403670,0.95242788,1.2319173,1.6674189,2.3527875,3.4499109,5.2424031,8.2360444,13.349647,22.283674,38.243824,67.384374,121.73344,225.21478,426.23651};
                T f1=(T)fep[this->m_dim-1];//Epanechikov
                cout<<"Normalization constant of Kernel type "<< this->TypeOfKernel<<": "<<f1<<endl;
                int i=0;
                for(i=0;i<=this->KERNEL_TABLE+1;i++)
                        {
                        this->KernelRad[i] = ((T)i)/(T)this->KERNEL_TABLE;
                        this->Kernel[i] = f1 *(1-this->KernelRad[i]*this->KernelRad[i]);
                        }
                this->Kernel[this->KERNEL_TABLE+1] =(T)0.0;
                };
        public:
        virtual T W(T u){
                int k;
                T wk=0.0;
                if(u<(T)1.0)
                        {
                        k = (int)(u*this->KERNEL_TABLE);
                        wk =( this->Kernel[k]  + (this->Kernel[k+1]-this->Kernel[k])*(u-this->KernelRad[k])*this->KERNEL_TABLE);
                        //rhoxyz+=wk*Part[pqStartA[i].p].Mass;
                        }
                return wk;
                };

        };




#endif

