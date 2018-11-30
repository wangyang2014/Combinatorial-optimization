/*
 * ==================================================================
 *         GLOBAL VARIABLES
 * ==================================================================
 */
string      instance;               //The name of the instance treated by the algorithm (no "/")
string      outputFolder;
ofstream    texFile;
long        timeMax;
int         srandSeed=0;



//MAIN VARIABLES ENCODING THE SOLUTION
int         n;                     //Number of vertices
int*        x;                     //The main vector of selected vertices
double**    m;                     //The matrix
int         k;                     //Number of points to be selected
int         h=0;                   //No vertex selected yet;


double      minDist    = INT_MAX,
            sumDist    = 0,
            bestKnown  = INT_MAX;
int         nminDist   = 0;


double*     minx;                       //The minimum distance from any vertex to a selected vertex
double*     sumx;                       //The sum for the same vertex
int*        nmin ;                      //How many times does it reach the minimum distance?

int*        iter;                       //the iteration numbers of the last change
long        it=1;                       //The current iteration, changed by function add
int*        dropTabu;
int*        addTabu;
int         genericDropTenure;
int         dropTenure;
int         genericAddTenure;
int         addTenure;
int         maxNoGain;                  //Maximum limit of iterations without any gain
int         noGain;                     //The current number of iterations without gain
int         maxIter;                    //iterMax*n = The maximum allowed number of iterations producing no new
                                        //better solution (compared to the best solution found so far in the search)
double      prevIterMinDist;
double      maxEverThisCycle=0;         //A cycle is performed in PRINT_RATE iterations
double      maxEverPrevCycle=0;
double      maxSumEverThisCycle=0;      //A cycle is performed in PRINT_RATE iterations
double      maxSumEverPrevCycle=0;

int         constantCycles=0;


double      epsilon=0.0001;             //used for some double comparisons
                                        //Think of <std::numeric_limits<double>::epsilon()
                                        //the smallest such that 1+epsilon!=1
                                        //attention: mmdp.cfg can change epsilon

clock_t     tStart;                     //The moment the algorithm started


std::map<string,string> configs; //configuration parameters


int recomputed;                         //The number of calls of linear case (3) of Algorithm 4 from the paper (the only one that might
                                        //take quadratic complexity if called many times per iteration


int * bestX=NULL;
double bestMinDist=-1;
double bestSumDist;
int lastIterWithNewBestPRecord = 0;        //The last iteration when an improvement was made
int timeAtNewBetter            = 0;        //Time passed at that moment, when new solution is found
int maxItersNoNewBetter        = 0;        //The maximum number of iterations with no new better solution encountered (so far)
