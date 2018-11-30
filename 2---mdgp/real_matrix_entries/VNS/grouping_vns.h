#include <stdio.h>

// Variable neighborhood search
#define DISTANCE_MIN                1
#define DISTANCE_MAX_COEF         0.3
#define DISTANCE_STEP               1

// Statistics
#define ST_COUNT                    8

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
#define bcon(X,Y) *((pinst+X)->bcon+Y)
#define bsol(X,Y) *((pgr+X)->bsol+Y)


// Structures
typedef struct
     {double *d;    /*  */
      double *con;  /*  */
      double *bcon; /*  */
     }Instance;

typedef struct
     {int s;        /*  */
      int bs;       /*  */
      int a;        /*  */
      int p;        /*  */
      int pg;       /*  */
      int sel;      /*  */
      int ob;       /*  */
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
     {double value;         /* solution value                      */
      double time_to_best;  /* time to the VNS best solution, secs */
      double total_time;    /* total time, secs                    */
     }Results;