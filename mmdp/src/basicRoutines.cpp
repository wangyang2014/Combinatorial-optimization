#include "basicRoutines.h"

#define freeMem(x) if(x) delete x;
void startTimer(){
    tStart=clock();
}
long timeCpuPassedMSec(){
	return (clock()-tStart)/((CLOCKS_PER_SEC/1000));
}
//*/

/*   //My checks tell that the above and below time computing procedures are equivalent.
     //On taurus I could never obtain a time of less than 10 milliseconds
#include <sys/times.h>
long clock2() {
  tms myTimes;
  long       ticks;
  double     msec;
  times( &myTimes );
  ticks = myTimes.tms_utime;
  msec = (double)ticks/(((double)sysconf(_SC_CLK_TCK))/1000.0);
  return (long)msec;
}

clock_t tStart;
void startTimer(){
    tStart=clock2();
}
long timeCpuPassedMSec(){
	return (clock2()-tStart);
}
//*/

double absolute(double x){if(x>0)return x;return -x;};

//0 = min, n-1 = max
int randomInt(int n){
     //double rnd01  ;
     return int(((double) n)*rand()/(RAND_MAX + 1.0));
     //rnd01  =  ((double)rand() )/ (RAND_MAX);
     //int toRet;
     //toRet = (int) floor(rnd01* n);
     //return  toRet;
}
//min = the minimum that can be generated, max-1 = the maximum that can be generated
int randomInt(int min, int max){
    return min+randomInt(max-min);
}
void switchInts(int& i1, int& i2){
	int tmp=i1;
	i1=i2;
	i2=tmp;
}


double keepTwoDecimals(double a){
    return 0.01*floor(0.5+a*100);
}
double keepFourDecimals(double a){
    return 0.0001*floor(0.5+a*10000);
}

/*  //Use these functions in case of various name conflicts (i.e. in Visual C++, header Windows.h)
template <class myType>
myType max(myType x, myType y){
	if(x>y)
		return x;
	return y;
}
template <class myType>
myType min(myType x, myType y){
	if(x<y)
		return x;
	return y;
}
//*/

//This routine returns 1 if and only if the string match over the last min(strlen(f*)) characters
//where f* is the smallest string
int sameWordEnding(char* f1,char*f2){
	int lenMin=strlen(f1);
	if(strlen(f2)<strlen(f1))
		lenMin = strlen(f2);
	f1 = f1 + strlen(f1)-lenMin;
	f2 = f2 + strlen(f2)-lenMin;
	return !strcmp(f1,f2);
}


