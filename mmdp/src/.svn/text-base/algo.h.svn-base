#ifndef ALGO_H_INCLUDED
#define ALGO_H_INCLUDED

#include "stdHeaders.h"
#include "basicRoutines.h"


#define CHECKRATE 100000

void setGenericTabuTenure(void);
void computeDropTenure(void);
void computeAddTenure(void);
void add(int vtx);
void drop(int vtx);
void recomputeMinDistAndNMin(int vtx);
void updateWhenDropped(int oldVtx);
void updateWhenAdded(int newVtx);
int acceptableDropMove(int chosenVtx);
void dropVertexOldest(void);
int acceptableAddMove(int chosenVtx, int &asspiredStatus);
int chooseAddGreedyMaxMin(void);
int firstPointFurthestFromAll(void);
void firstTwo(void);
void addVertexConstructive(void);
void constructiveAlg(int level);
void constructiveAlg(void);
void dropAddLocalSearch(void);
void dropAddAlg(int level);
void dropAddAlg(void);
int checkMMDP(int*,double);
int checkMDP(int*,double&);

int checkNewBetterSolution();
double bestMinDistUntilNow();
int itersNoNewBetter();







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





#endif // ALGO_H_INCLUDED
