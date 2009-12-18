//last change Arman Khalatyan akhalatyan@aip.de
//2004-06-18
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

//int __gxx_personality_v0;

//    void    starttimer_(void);
//    float  stoptimer_(int  *flag);
//    float  gettimer_(void);
//    double drand_(int);


clock_t startm, stopm;
void  starttimer_( void)
    //  Start Timer
{
    
    if ( (startm = clock()) == -1) {printf("Error calling clock\n\n");exit(1);}
}
float   stoptimer_( int *flag)
    //Stop and print
{
    //printf("flag = %d", *flag);
    if ( (stopm = clock()) == -1) {printf("Error calling clock\n\n");exit(1);}
    
    float run_time = ((float)(stopm-startm))/CLOCKS_PER_SEC/60.0f;
    if(*flag!=0)
      printf( "\n%6.4f min used by the processor.\n", run_time);
   
    return run_time;
}

float   gettimer_( void)
    //Stop and print
{
       
    float timer;
    if ( (timer = float(clock())) == -1) 
    {
	printf("Error calling clock\n\n");
	exit(1);
    }

    // printf( "\nclock = %6.4f\n", timer);
   
    return timer/CLOCKS_PER_SEC;
}
double drand_(int seed)
{
 srand(seed);
 return double(rand())/float(RAND_MAX);
}

float mydrand48(void)
{
float x = rand()/float(RAND_MAX);
return x;

};

