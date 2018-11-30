#include <stdio.h>

// Genetic algorithm
#define POPULATION_SIZE           100

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
#define pop(X,Y) *((psol+Y)->pop+X)
#define pfv(X) *(psol->pfv+X)
#define bsol(X,Y) *((pgr+X)->bsol+Y)
#define popsz(X,Y) *((pgr+X)->popsz+Y)
#define oob(X,Y) *((pgr+X)->oob+Y)
#define ob1(X,Y) *((pgr+X)->ob1+Y)
#define ob2(X,Y) *((pgr+X)->ob2+Y)
#define is(X,Y) *((pgr+X)->is+Y)
#define p(X,Y) *((pgr+X)->p+Y)


// Structures
typedef struct
     {double *d;    /*  */
      double *con;  /*  */
     }Instance;

typedef struct
     {int *pop;     /*  */
      double *pfv;  /*  */
      int s;        /*  */
      int bs;       /*  */
      int a;        /*  */
      int c;        /*  */
      int p;        /*  */
      int pg;       /*  */
      int m;        /*  */
      int sel;      /*  */
      int r;        /*  */
      long perf;    /*  */
     }Solution;

typedef struct
     {int *bsol;    /*  */
      int *popsz;   /*  */
      int *oob;     /*  */
      int *ob1;     /*  */
      int *ob2;     /*  */
      int *is;      /*  */
      int *p;       /*  */
      int sz;       /*  */
      int bsz;      /*  */
      int mins;     /*  */
      int maxs;     /*  */
      int oc;       /*  */
      int c1;       /*  */
      int c2;       /*  */
      int i1;       /*  */
      int i2;       /*  */
      int cand;     /*  */
      int c;        /*  */
     }Group;

typedef struct
     {double value;         /* solution value                      */
      double time_to_best;  /* time to the GA best solution, secs  */
      double total_time;    /* total time, secs                    */
     }Results;