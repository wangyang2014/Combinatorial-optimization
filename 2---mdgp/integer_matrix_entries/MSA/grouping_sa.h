#include <stdio.h>

// Simulated annealing
#define COOLING_FACTOR           0.95
#define MIN_TEMPERATURE        0.0001
#define REPETITION_COEF           300
#define SAMPLE_SIZE              1000
#define RELOCATING_PROBAB         0.3

// Statistics
#define ST_COUNT                   10

// Memory allocation
#define ALS(X,Y,Z) if ((X=(Y *)calloc(Z,sizeof(Y)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALI(X,Z) if ((X=(int *)calloc(Z,sizeof(int)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}

// Matrices
#define d(X,Y) *((pinst+X)->d+Y)
#define con(X,Y) *((pinst+X)->con+Y)
#define bsol(X,Y) *((pgr+X)->bsol+Y)


// Structures
typedef struct
     {int *d;       /*  */
      int *con;     /*  */
     }Instance;

typedef struct
     {int s;        /*  */
      int bs;       /*  */
      int a;        /*  */
      int p;        /*  */
      long perf;    /*  */
     }Solution;

typedef struct
     {int sz;       /*  */
      int bsz;      /*  */
      int mins;     /*  */
      int maxs;     /*  */
      int c;        /*  */
      int *bsol;    /*  */
     }Group;

typedef struct
     {long value;           /* solution value                      */
      double time_to_best;  /* time to the SA best solution, secs  */
      double total_time;    /* total time, secs                    */
     }Results;