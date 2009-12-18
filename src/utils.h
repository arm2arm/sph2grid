#ifndef _CRange_
#define _CRange_
#include <iostream>
#include <algorithm>
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
void print(char *str){
	cout<<str<<" Min="<<Min<<"\tMax="<<Max<<"\t mean="<<m_sum/m_np<<endl;

};
};

#define printOpenGLError() printOglError(__FILE__, __LINE__)

int printOglError(char *file, int line);
float   stoptimer_( int *flag);
void  starttimer_( void);
float mydrand48(void);
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795f
#endif

#endif

