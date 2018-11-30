#ifndef BESTSOLANDTEXOUTPUT_H_INCLUDED
#define BESTSOLANDTEXOUTPUT_H_INCLUDED

#include "stdHeaders.h"
#include "basicRoutines.h"


void printMatrix(ostream&,int* );
int itersNoNewBetter(void);
int checkNewBetterSolution(void);
void nextLevel(void);
void resetBestSolution(void);
int checkMMDP(int *x, double minDist);
int checkMDP(int *x, double &sumDist);
void saveThisSpecificSolution(int *x, double minDist, int keepOutputFolderInFilename);
void saveThisSpecificSolution(int *x, double minDist);
void saveBestSolutionBetterPerf(void);
void saveBestSolutionEqualPerf(void);
void finishCheckAndWriteTex(int *x, double minDist, double sumDist);
void finishCheckAndWriteTex(void);
void printBestSolution(void);
double bestMinDistUntilNow(void);


extern int *  bestX ;
extern double bestMinDist;
extern double bestSumDist;
extern int lastIterWithNewBestPRecord ;        //The last iteration when an improvement was made
extern int timeAtNewBetter            ;        //Time passed at that moment, when new solution is found
extern int maxItersNoNewBetter        ;        //The maximum number of iterations with no new better solution encountered (so far)


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


extern int         genericDropTenure;
extern int         dropTenure;
extern int         genericAddTenure;
extern int         addTenure;
extern int         h;                         //No vertex selected yet;

extern int         recomputed;                //The number of calls of linear case (3) of Algorithm 4 from the paper (the only one that might
                                              //take quadratic complexity if called many times per iteration

extern double      minDist    ;
extern double      sumDist    ;
extern double      bestKnown  ;
extern int         nminDist   ;



#endif // BESTSOLANDTEXOUTPUT_H_INCLUDED
