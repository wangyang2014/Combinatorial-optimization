/* Program: an implementation of the variable neighborhood search algorithm 
   for the maximally diverse grouping problem. 
   Author: Gintaras Palubeckis
   Date: 2011-07-29
   Language: C
   Some of the input data are supplied through the parameters and the rest 
   through the (input) file. The program terminates when a specified time 
   limit is reached. 
   Parameters:
     - input file name;
     - output file name;
     - seed for random number generator;
     - time limit (in seconds);
     - a pointer to structure 'Results' for writing down the performance 
       characteristics. The structure is defined in file 'grouping_vns.h'. 
       It is possible to have the last parameter null: in this case, 
       information is directed only to the output file.
   Examples of invocation:
       either
     Results *pres;
     char in_file_name[80],out_file_name[80];
     double seed=1000;
     pres=(Results *)calloc(1,sizeof(Results));
     strcpy(in_file_name,"D:\\data\\RanReal_n012_ds_06.txt");
     strcpy(out_file_name,"D:\\temp\\RanReal_n012_ds_06.res");
     grouping_vns(in_file_name,out_file_name,seed,20,pres);
       or just
     char in_file_name[80],out_file_name[80];
     double seed=1000;
     strcpy(in_file_name,"D:\\data\\RanReal_n012_ds_06.txt");
     strcpy(out_file_name,"D:\\temp\\RanReal_n012_ds_06.res");
     grouping_vns(in_file_name,out_file_name,seed,20,NULL);
   Input file contains:
     - the number of elements n;
     - the number of groups gc;
     - the string constant either "ss" or "ds"; in the former case, 
       all groups are of the same size, whereas in the second case, the size of 
       groups is allowed to vary within specified bounds;
     - for each group, an integer specifying the minimum possible size followed 
       by an integer specifying the maximum possible size; 
     - for each pair of elements, the triplet in which the first two integers 
       denote elements themselves and the third number expresses the dissimilarity 
       between them (the elements in the file are numbered starting from 0).
   Example of the input file:
   12 4 ds 3 4 3 3 2 4 3 4 
   0 1 60.709
   0 2 99.446
   0 3 25.333
   0 4 91.585
   0 5 69.584
   0 6 6.232
   0 7 95.986
   0 8 44.604
   0 9 37.539
   0 10 70.317
   0 11 36.884
   1 2 94.406
   .............
   .............
   .............
*/



#include <alloc.h>
#include <process.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "grouping_vns.h"



 double random(double *seed,double coef)
 {double rd,rf;
  rd=16807*(*seed);rf=floor(rd/coef);
  *seed=rd-rf*coef;
  return(*seed/(coef+1));
 }


 double take_time(int *time_values,clock_t start)
 {int i;
  int hours,mins;
  long longsecs;
  double elapsed_sec;
  clock_t end;
  end=clock();
  elapsed_sec=(end-start)/CLK_TCK;
  longsecs=elapsed_sec;
  for (i=1;i<=4;i++) time_values[i]=0;
  hours=(int)(longsecs/3600);
  if (hours>0)     /* more than an hour  */
     {time_values[1]=hours;longsecs-=hours*3600;}
  mins=(int)(longsecs/60);
  if (mins>0)     /* more than a minute */
     {time_values[2]=mins;longsecs-=mins*60;}
  time_values[3]=(int)longsecs;
  time_values[4]=elapsed_sec*1000-(long)elapsed_sec*1000;
  return elapsed_sec;
 }


 void check_best_solution(FILE *out,int n,int gc,double best_value,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i,j,k;
  int count;
  double value_obj=0,value_groups=0;

  for (k=1;k<=gc;k++) if ((pgr+k)->bsz<(pgr+k)->mins || (pgr+k)->bsz>(pgr+k)->maxs)
     {fprintf(out,
          "!!! group size violation in VNS solution: %4d   %4d\n",
          k,(pgr+k)->bsz);
      exit(1);
     }
  for (k=1;k<=gc;k++)
     {count=0;
      for (i=1;i<=n;i++) if ((psol+i)->bs==k) count++;
      if (count<(pgr+k)->mins || count>(pgr+k)->maxs)
         {fprintf(out,
          "!!! group size violation in VNS solution (2): %4d   %4d\n",
              k,count);
          exit(1);
         }
     }
  for (i=1;i<n;i++) for (j=i+1;j<=n;j++)
     {if ((psol+i)->bs!=(psol+j)->bs) continue;
      value_obj+=(d(i,j));
     }
  if (value_obj<best_value-0.00001 || value_obj>best_value+0.00001)
     {fprintf(out,
          "!!! some discrepancy in VNS solution values (1): %10lf   %10lf\n",
          best_value,value_obj);
      exit(1);
     }
  for (k=1;k<=gc;k++)
  for (i=1;i<(pgr+k)->bsz;i++) for (j=i+1;j<=(pgr+k)->bsz;j++)
      value_groups+=(d(bsol(k,i),bsol(k,j)));
  if (value_groups<best_value-0.00001 || value_groups>best_value+0.00001)
     {fprintf(out,
          "!!! some discrepancy in VNS solution values (2): %10lf   %10lf\n",
          best_value,value_groups);
      exit(1);
     }
 }


 double random_start(int n,int gc,double coef,double *seed,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i,j,k,m;
  int r=0;
  int kg,mg;
  double sol_value=0;

// Randomly generating permutation
  for (i=1;i<=n;i++)
     {(psol+i)->a=i;
      for (j=1;j<=gc;j++) con(i,j)=0;
     }
  for (i=1;i<=n;i++)
     {k=random(seed,coef)*(n-i+1);k+=i;
      (psol+i)->p=(psol+k)->a;(psol+k)->a=(psol+i)->a;
     }

// Constructing solution
  for (j=1;j<=gc;j++)
     {for (i=1;i<=(pgr+j)->mins;i++)
         {r++;
          k=(psol+r)->p;
          (psol+k)->s=(psol+k)->bs=j;
         }
      (pgr+j)->sz=(pgr+j)->mins;
     }
  if (r<n) for (j=1;j<=gc;j++)
     {for (i=(pgr+j)->mins+1;i<=(pgr+j)->maxs;i++)
         {r++;
          k=(psol+r)->p;
          (psol+k)->s=(psol+k)->bs=j;((pgr+j)->sz)++;
          if (r>=n) break;
         }
      if (r>=n) break;
     }
// Calculating the value of the objective function
  for (k=1;k<n;k++)
     {kg=(psol+k)->s;
      for (m=k+1;m<=n;m++)
         {mg=(psol+m)->s;
          (con(k,mg))+=(d(k,m));(con(m,kg))+=(d(k,m));
          if (kg==mg) sol_value+=(d(k,m));
         }
     }
  return sol_value;
 }


 void update_temp_solution(int n,int k,int m,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i;
  int kg,mg;

  if (m>0)
     {kg=(psol+k)->s;mg=(psol+m)->s;
      for (i=1;i<=n;i++)
         {if (i==k || i==m) continue;
          con(i,kg)+=(d(i,m)-d(i,k));
          con(i,mg)+=(d(i,k)-d(i,m));
         }
      con(k,kg)+=(d(k,m));con(k,mg)-=(d(k,m));
      con(m,kg)-=(d(k,m));con(m,mg)+=(d(k,m));
      (psol+k)->s=mg;(psol+m)->s=kg;
     }
   else
     {kg=(psol+k)->s;mg=-m;
      for (i=1;i<=n;i++)
         {if (i==k) continue;
          con(i,kg)-=(d(i,k));
          con(i,mg)+=(d(i,k));
         }
      (psol+k)->s=mg;
      ((pgr+kg)->sz)--;((pgr+mg)->sz)++;
     }
 }


 double shake(int n,int gc,int dist,double sol_value,double coef,
      double *seed,Instance *pinst,Solution *psol,Group *pgr)
 {int i,j,k,l,m;
  int kg,mg;
  int swap_count=0,swap_req;
  int reloc_count=0,reloc_req;
  int ind,ind1,ind2;
  int rem,count;
  double dif;

  for (i=1;i<=n;i++)
     {(psol+i)->s=(psol+i)->bs;(psol+i)->sel=0;(psol+i)->ob=i;
      for (j=1;j<=gc;j++) con(i,j)=bcon(i,j);
     }
  for (j=1;j<=gc;j++) (pgr+j)->sz=(pgr+j)->bsz;
  rem=n;
  swap_req=dist/2;

  while (swap_count<swap_req)
     {ind1=random(seed,coef)*rem+1;
      k=(psol+ind1)->ob;kg=(psol+k)->s;
      for (i=1,count=0;i<=rem;i++)
          if (i!=ind1 && (psol+(psol+i)->ob)->s!=kg)
             {count++;(psol+count)->p=i;}
      if (count<=0) break;
      ind2=random(seed,coef)*count+1;ind2=(psol+ind2)->p;
      m=(psol+ind2)->ob;mg=(psol+m)->s;
      swap_count++;
      dif=con(k,mg)-con(k,kg)+con(m,kg)-con(m,mg)-2*d(k,m);
      sol_value+=dif;
      update_temp_solution(n,k,m,pinst,psol,pgr);
      (psol+k)->sel=2;(psol+m)->sel=2;
      if (swap_count>=swap_req || rem<4) break;
      if (ind2<rem)
         {(psol+ind1)->ob=(psol+rem)->ob;rem--;
          (psol+ind2)->ob=(psol+rem)->ob;rem--;
         }
       else
         {rem--;(psol+ind1)->ob=(psol+rem)->ob;rem--;}
     }
  reloc_req=dist-2*swap_count;
  if (reloc_req<=0) return sol_value;
  for (i=1,rem=0;i<=n;i++)
     {if ((psol+i)->sel>1 || (pgr+(psol+i)->s)->sz<=(pgr+(psol+i)->s)->mins) continue;
      rem++;(psol+rem)->ob=i;
     }
  if (rem<=0) return sol_value;
  while (reloc_count<reloc_req)
     {ind=random(seed,coef)*rem+1;
      k=(psol+ind)->ob;kg=(psol+k)->s;
      for (i=1,count=0;i<=gc;i++)
         {if (i==kg || (pgr+i)->sz>=(pgr+i)->maxs) continue;
          count++;(psol+count)->p=i;
         }
      if (count==0) break;
      l=random(seed,coef)*count+1;j=(psol+l)->p;
      dif=con(k,j)-con(k,kg);
      sol_value+=dif;
      update_temp_solution(n,k,-j,pinst,psol,pgr);
      reloc_count++;
      if (reloc_count>=reloc_req) break;
      for (i=1,count=0;i<=rem;i++)
         {m=(psol+i)->ob;
          if (i==ind || (pgr+(psol+m)->s)->sz<=(pgr+(psol+m)->s)->mins) continue;
          count++;(psol+count)->ob=m;
         }
      if (count<=0) break;
      rem=count;
     }
  return sol_value;
 }


 double local_search(int n,int gc,double best_value,double coef,
      double *seed,Instance *pinst,Solution *psol,Group *pgr)
   {int i,j,j0,k,k0,m,m0;
  int kg,mg;
  int ind1,ind2;
  double dif;
  double temp_best=1;
  double sol_value;

  sol_value=best_value;
// Randomly generating permutations
  for (i=1;i<=n;i++) (psol+i)->a=i;
  for (i=1;i<=n;i++)
     {k=random(seed,coef)*(n-i+1);k+=i;
      (psol+i)->p=(psol+k)->a;(psol+k)->a=(psol+i)->a;
     }
  for (i=1;i<=gc;i++) (psol+i)->a=i;
  for (i=1;i<=gc;i++)
     {k=random(seed,coef)*(gc-i+1);k+=i;
      (psol+i)->pg=(psol+k)->a;(psol+k)->a=(psol+i)->a;
     }

  while (temp_best>0)
     {temp_best=0;
      for (k0=1;k0<=n;k0++)
         {k=(psol+k0)->p;kg=(psol+k)->s;
          if ((pgr+kg)->sz>(pgr+kg)->mins)
// Trying to relocate the element under consideration
             {for (j0=1;j0<=gc;j0++)
                 {j=(psol+j0)->pg;if (j==kg || (pgr+j)->sz>=(pgr+j)->maxs) continue;
                  dif=con(k,j)-con(k,kg);
                  if (dif>0.00001) {temp_best=dif;ind1=k;ind2=-j;break;}
                 } // for j0
              if (temp_best>0.00001) break;
             }
// All pairs of elements in different groups are searched. 
// Looking for the best pairwise exchange of elements 
// (these elements are ind1, ind2)
          for (m0=k0+1;m0<=n;m0++)
             {m=(psol+m0)->p;if ((psol+m)->s==kg) continue;
              mg=(psol+m)->s;
              dif=con(k,mg)-con(k,kg)+con(m,kg)-con(m,mg)-2*d(k,m);
// If "dif>0", then better solution has been found
              if (dif>0.00001) {temp_best=dif;ind1=k;ind2=m;break;}
             } // for m0
          if (temp_best>0.00001) break;
         } // for k0
      if (temp_best>0.00001)
// Improving move is detected
         {update_temp_solution(n,ind1,ind2,pinst,psol,pgr);
// Solution value is updated
          sol_value+=temp_best;
         }
     } // while
  return sol_value;
 }


 void store_new_best_sol(int n,int gc,int st,int *time_values_opt,
      double *time_to_opt,clock_t start_time,Instance *pinst,
      Solution *psol,Group *pgr)
 {int i,j;

  for (j=1;j<=gc;j++) {(pgr+j)->bsz=(pgr+j)->sz;(pgr+j)->c=0;}
  for (i=1;i<=n;i++)
     {(psol+i)->bs=(psol+i)->s;
      j=(psol+i)->s;((pgr+j)->c)++;bsol(j,(pgr+j)->c)=i;
      for (j=1;j<=gc;j++) bcon(i,j)=con(i,j);
     }
  ((psol+1)->perf)++;(psol+2)->perf=st;
  *time_to_opt=take_time(time_values_opt,start_time);
 }


 double variable_neighborhood_search(int n,int gc,int dist_min,
      int dist_max,int dist_step,long time_limit,
      int *time_values_opt,double *time_to_opt,double *seed1,
      clock_t start_time,Instance *pinst,Solution *psol,Group *pgr)
 {int dist;
  int st=0;
  int stop_cond=0;
  double best_value,sol_value;
  double seed2,seed3,coef;
  double elapsed_time;
  clock_t end;

  coef=2048;coef*=1024;coef*=1024;coef-=1;
  seed2=2*(*seed1);seed3=3*(*seed1);
  (psol+1)->perf=0;

// Generating initial solution randomly
  sol_value=random_start(n,gc,coef,seed1,pinst,psol,pgr);
  sol_value=local_search(n,gc,sol_value,coef,&seed3,pinst,psol,pgr);
  store_new_best_sol(n,gc,0,time_values_opt,time_to_opt,start_time,
      pinst,psol,pgr);
  best_value=sol_value;(psol+1)->perf=0;
  while (stop_cond==0)
     {dist=dist_min;
      st++;
      while (dist<=dist_max && stop_cond==0)
         {sol_value=shake(n,gc,dist,best_value,coef,&seed2,pinst,psol,pgr);
          sol_value=local_search(n,gc,sol_value,coef,&seed3,pinst,psol,pgr);
          if (sol_value>best_value)
             {store_new_best_sol(n,gc,st,time_values_opt,
                  time_to_opt,start_time,pinst,psol,pgr);
              best_value=sol_value;
              dist=dist_min;
             }
           else dist+=dist_step;
// Checking variable neighborhood search termination rule
          end=clock();
          elapsed_time=(end-start_time)/CLK_TCK;
          if (elapsed_time>=time_limit) stop_cond=1;
         } // inner while
     } // outer while
// Storing the total number of variable neighborhood search starts executed
  (psol+3)->perf=st;
  return best_value;
 }





// 'grouping' is the VNS procedure for creating maximally diverse groups from a given set of elements
 void grouping_vns(char *in_file_name,char *out_file_name,double seed,
      long time_limit,Results *pres)
 {FILE *out,*in;
  Instance *pinst;
  Solution *psol;
  Group *pgr;

  int i,j,k,p,r;
  int n,gc;
  int dist_max;
  int minsize,maxsize=-1;
  int lo;
  int time_values[5],time_values_opt[5];
  double best_value;
  double we;
  double time_in_seconds,time_to_opt;
  double seed_saved;
  clock_t start_time;
  char str[5];

  if ((in=fopen(in_file_name,"r"))==NULL)
     {printf("  fopen failed for input");exit(1);}
// n is the number of elements
  fscanf(in,"%d %d %s",&n,&gc,&str);
// Alternative input
  //fscanf(in,"%d",&n);
  if ((out=fopen(out_file_name,"w"))==NULL)
     {printf("  fopen failed for output  %s",out_file_name);exit(1);}
  seed_saved=seed;
// Allocation of core memory for an array of Instance structure objects 
// that contains the dissimilarity matrix
  ALS(pinst,Instance,n+1)
  for (i=0;i<=n;i++) ALD((pinst+i)->d,n+1)
  for (i=0;i<=n;i++) ALD((pinst+i)->con,gc+1)
  for (i=0;i<=n;i++) ALD((pinst+i)->bcon,gc+1)
// Allocation of core memory for an array of Solution structure objects 
// that contains data used in various methods. 
// ST_COUNT is the length of performance statistics array 'perf'
  if (n>ST_COUNT) i=n; else i=ST_COUNT;
  ALS(psol,Solution,i+1)
// Allocation of core memory for an array of Group structure objects 
// that contains data used in the construction of groups
  ALS(pgr,Group,gc+1)
  minsize=n+1;
  for (i=1;i<=gc;i++)
     {fscanf(in,"%d %d",&p,&r);
      (pgr+i)->mins=p;(pgr+i)->maxs=r;
      if (r>maxsize) maxsize=r;if (p<minsize) minsize=p;
     }
  for (i=0;i<=gc;i++) ALI((pgr+i)->bsol,maxsize+1)

// Input of the dissimilarity matrix
  lo=(n*(n-1))/2;
  for (k=1;k<=lo;k++)
     {fscanf(in,"%d %d %lf",&i,&j,&we);i++;j++;
      d(i,j)=we;d(j,i)=we;
     }
  for (i=1;i<=n;i++) d(i,i)=0;
// Alternative input
/*  for (i=1;i<=n;i++) for (j=1;j<=n;j++)
     {fscanf(in,"%lf",&we);
      if (i==j) d(i,j)=0; else d(i,j)=we;
     }*/

// Preparation of the VNS parameter "dist_max"
  dist_max=n*DISTANCE_MAX_COEF;

  start_time=clock();
// Running variable neighborhood search
  best_value=variable_neighborhood_search(n,gc,DISTANCE_MIN,
      dist_max,DISTANCE_STEP,time_limit,time_values_opt,
      &time_to_opt,&seed,start_time,pinst,psol,pgr);
  time_in_seconds=take_time(time_values,start_time);
  check_best_solution(out,n,gc,best_value,pinst,psol,pgr);

// Printing parameters, data characteristics and some statistics 
// regarding the performance of the algorithm
  fprintf(out,"   parameters:                                 \n");
  fprintf(out,"      DISTANCE_MIN                           = %5d\n",
      DISTANCE_MIN);
  fprintf(out,"      DISTANCE_MAX_COEF                      = %6lf\n",
      DISTANCE_MAX_COEF);
  fprintf(out,"      dist_max                               = %5d\n",
      dist_max);
  fprintf(out,"      DISTANCE_STEP                          = %5d\n",
      DISTANCE_STEP);
  fprintf(out,"      seed for random number generator       = %8lf\n",
      seed_saved);
  fprintf(out,"   number of elements                         = %5d\n",n);
  fprintf(out,"   number of groups                          = %5d\n",gc);
  fprintf(out,"   minimum group size                        = %5d\n",minsize);
  fprintf(out,"   maximum group size                        = %5d\n",maxsize);
  fprintf(out,"   time limit                                = %5ld\n",
      time_limit);
  fprintf(out,"   number of VNS starts executed             = %4ld\n",
      (psol+3)->perf);
  fprintf(out,"   number of improvements during VNS         = %4ld\n",
      (psol+1)->perf);
  fprintf(out,"   last improvement at VNS start no.         = %3ld\n",
      (psol+2)->perf);
  fprintf(out,"   value found by VNS                        = %8lf\n",
      best_value);
  (psol+4)->perf=(long)time_to_opt;
  fprintf(out,"   time to the best solution: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values_opt[1],time_values_opt[2],time_values_opt[3],
      time_values_opt[4],(long)time_to_opt);
  (psol+5)->perf=(long)time_in_seconds;
  fprintf(out,"   total time of VNS: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values[1],time_values[2],time_values[3],
      time_values[4],(long)time_in_seconds);
  fflush(out);

  if (pres!=NULL)
     {pres->value=best_value;
      pres->total_time=time_in_seconds;
      pres->time_to_best=time_to_opt;
     }

// Releasing the memory
  for (i=0;i<=gc;i++) if ((pgr+i)->bsol!=NULL) free((pgr+i)->bsol);
  if (pgr!=NULL) free(pgr);
  if (psol!=NULL) free(psol);
  for (i=0;i<=n;i++) if ((pinst+i)->bcon!=NULL) free((pinst+i)->bcon);
  for (i=0;i<=n;i++) if ((pinst+i)->con!=NULL) free((pinst+i)->con);
  for (i=0;i<=n;i++) if ((pinst+i)->d!=NULL) free((pinst+i)->d);
  if (pinst!=NULL) free(pinst);
// Closing input and output files
  fclose(out);fclose(in);
 }
