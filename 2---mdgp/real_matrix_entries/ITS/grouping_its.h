#include <stdio.h>

// Iterated tabu search
#define TABU_TIME                  10
#define TABU_COEF                   4
#define PER_COEF1                 0.1
#define PER_COEF2                 0.5
#define MIN_PER_COUNT              10
#define LIST_SIZE1                 10//25
#define LIST_SIZE2                300
#define BUCKET_SIZE                10
#define REL_PROB                  0.4
#define PROBL_SIZE_CONST          300
#define ITER_COUNT_SMALLER        100
#define ITER_COUNT_BIGGER         200

// Statistics
#define ST_COUNT                    8

// Constants
#define POS_LARGE               30000
#define NEG_LARGE              -30000

// Memory allocation
#define ALS(X,Y,Z) if ((X=(Y *)calloc(Z,sizeof(Y)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALI(X,Z) if ((X=(int *)calloc(Z,sizeof(int)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALD(X,Z) if ((X=(double *)calloc(Z,sizeof(double)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}

// Matrices
#define d(X,Y) *((pinst+X)->d+Y)
#define con(X,Y) *((pinst+X)->con+Y)
#define texch(X,Y) *((psol+X)->texch+Y)
#define trel(X,Y) *((psol+X)->trel+Y)
#define bsol(X,Y) *((pgr+X)->bsol+Y)


// Structures
typedef struct
     {double *d;    /*  */
      double *con;  /*  */
     }Instance;

typedef struct
     {int *texch;   /*  */
      int *trel;    /*  */
      int s;        /*  */
      int bs;       /*  */
      int a;        /*  */
      int p;        /*  */
      int pg;       /*  */
      int sel;      /*  */
      int t1;       /*  */
      int t2;       /*  */
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
     {double minv;  /*  */
      double cr;    /*  */
      int v1;       /*  */
      int v2;       /*  */
      int bl;       /*  */
      int bu;       /*  */
      int ind;      /*  */
     }Move;

typedef struct
     {double value;         /* solution value                      */
      double time_to_best;  /* time to the ITS best solution, secs */
      double total_time;    /* total time, secs                    */
     }Results;