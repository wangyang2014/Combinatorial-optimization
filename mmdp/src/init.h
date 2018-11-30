#ifndef INIT_H_INCLUDED
#define INIT_H_INCLUDED
#include "stdHeaders.h"


void init(int , char**, string );
void usage();
void readGlover(ifstream& ,int);
void readRun(ifstream& );
void readMatrix(char* );
void getInstanceNameAndK(char*,string&);
void readConfigFile(const char*);
void readConfigParams();


extern std::map<string,string> configs; //configuration parameters

extern string      instance;               //The name of the instance treated by the algorithm (no "/")
extern string      outputFolder;
extern ofstream    texFile;
extern long        timeMax;


extern int         n;                     //Number of vertices
extern int*        x;                     //The main vector of selected vertices
extern int         k;                     //Number of points to be selected
extern double**    m;                     //The matrix
extern int*        iter;                  //the iteration numbers of the last change
extern long        it;                       //The current iteration, changed by function add
extern double*     minx;                       //The minimum distance from any vertex to a selected vertex
extern double*     sumx;                       //The sum for the same vertex
extern int*        nmin ;                      //How many times does it reach the minimum distance?
extern int*        dropTabu;
extern int*        addTabu;
extern double      bestKnown  ;
extern double      epsilon ;              //used for some double comparisons
extern int         maxIter;                    //iterMax*n = The maximum allowed number of iterations producing no new
extern int         srandSeed;
extern int         maxNoGain;                  //Maximum limit of iterations without any gain


#endif // INIT_H_INCLUDED
