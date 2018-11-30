/* Program: an implementation of the simulated annealing algorithm 
   for the maximally diverse grouping problem. 
   Author: Gintaras Palubeckis
   Date: 2011-08-15
   Language: C
   Some of the input data are supplied through the parameters and the rest 
   through the (input) file. The program terminates either when a specified time 
   limit is reached or after completing the first run of the simulated 
   annealing (SA) procedure when the allotted time interval is insufficient 
   to do this. In the former case, SA procedure is restarted with a new 
   randomly generated solution.
   Parameters:
     - input file name;
     - output file name;
     - seed for random number generator;
     - time limit (in seconds);
     - a pointer to structure 'Results' for writing down the performance 
       characteristics. The structure is defined in file 'grouping_sa.h'. 
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
     grouping_sa(in_file_name,out_file_name,seed,20,pres);
       or just
     char in_file_name[80],out_file_name[80];
     double seed=1000;
     strcpy(in_file_name,"D:\\data\\RanReal_n012_ds_06.txt");
     strcpy(out_file_name,"D:\\temp\\RanReal_n012_ds_06.res");
     grouping_sa(in_file_name,out_file_name,seed,20,NULL);
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
#include "grouping_sa.h"



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
          "!!! group size violation in SA solution: %4d   %4d\n",
          k,(pgr+k)->bsz);
      exit(1);
     }
  for (k=1;k<=gc;k++)
     {count=0;
      for (i=1;i<=n;i++) if ((psol+i)->bs==k) count++;
      if (count<(pgr+k)->mins || count>(pgr+k)->maxs)
         {fprintf(out,
          "!!! group size violation in SA solution (2): %4d   %4d\n",
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
          "!!! some discrepancy in SA solution values (1): %10lf   %10lf\n",
          best_value,value_obj);
      exit(1);
     }
  for (k=1;k<=gc;k++)
  for (i=1;i<(pgr+k)->bsz;i++) for (j=i+1;j<=(pgr+k)->bsz;j++)
      value_groups+=(d(bsol(k,i),bsol(k,j)));
  if (value_groups<best_value-0.00001 || value_groups>best_value+0.00001)
     {fprintf(out,
          "!!! some discrepancy in SA solution values (2): %10lf   %10lf\n",
          best_value,value_groups);
      exit(1);
     }
 }


 double start_sol(int n,int gc,double coef,double *seed,
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
          k=(psol+r)->p;(psol+k)->s=j;
         }
      (pgr+j)->sz=(pgr+j)->mins;
     }
  if (r<n) for (j=1;j<=gc;j++)
     {for (i=(pgr+j)->mins+1;i<=(pgr+j)->maxs;i++)
         {r++;
          k=(psol+r)->p;
          (psol+k)->s=j;((pgr+j)->sz)++;
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


 void store_new_best_sol(int n,int gc,int st,int *time_values_opt,
      double *time_to_opt,clock_t start_time,Solution *psol,Group *pgr)
 {int i,j;

  for (j=1;j<=gc;j++) {(pgr+j)->bsz=(pgr+j)->sz;(pgr+j)->c=0;}
  for (i=1;i<=n;i++)
     {(psol+i)->bs=(psol+i)->s;
      j=(psol+i)->s;((pgr+j)->c)++;bsol(j,(pgr+j)->c)=i;
     }
  ((psol+1)->perf)++;(psol+2)->perf=st;
  *time_to_opt=take_time(time_values_opt,start_time);
 }


 double set_max_temperature(int n,int pair_exch_count,
      double coef,double *seed,Instance *pinst,Solution *psol)
 {int i,k,m;
  int kg,mg;
  double max_temp=-1,dif;

  for (i=1;i<=pair_exch_count;i++)
     {k=random(seed,coef)*n+1;kg=(psol+k)->s;mg=kg;
      while (mg==kg)
         {m=random(seed,coef)*n+1;
          mg=(psol+m)->s;
         }
      dif=con(k,mg)-con(k,kg)+con(m,kg)-con(m,mg)-2*d(k,m);
      if (dif<0) dif=-dif;
      if (dif>max_temp) max_temp=dif;
     }
  return max_temp;
 }


 double simulated_annealing(int n,int gc,int t_val_count,int rep_count,
      double best_value,double max_temp,long start,long time_limit,
      double coolf,double rel_prob,double coef,int *stop_cond,
      int *time_values_opt,double *time_to_opt,double *seed,
      clock_t start_time,Instance *pinst,Solution *psol,Group *pgr)
 {int j,k,m,r,s;
  int kg,mg;
  int size;
  int reloc;
  long elapsed_time;
  double sol_value;
  double dif;
  double t;
  double acc_val,db;
  clock_t end_time;

  if (start>1) sol_value=start_sol(n,gc,coef,seed,pinst,psol,pgr);
   else sol_value=best_value;
  t=max_temp;

  for (r=1;r<=t_val_count;r++)
     {for (s=1;s<=rep_count;s++)
         {if (rel_prob>0)
             {db=random(seed,coef);
              if (db<=rel_prob) reloc=1; else reloc=0;
             }
           else reloc=0;
          if (reloc==0)
             {k=random(seed,coef)*n+1;kg=(psol+k)->s;mg=kg;
              while (mg==kg)
                 {m=random(seed,coef)*n+1;
                  mg=(psol+m)->s;
                 }
              dif=con(k,mg)-con(k,kg)+con(m,kg)-con(m,mg)-2*d(k,m);
              if (dif>=0)
                 {update_temp_solution(n,k,m,pinst,psol,pgr);
                  sol_value+=dif;
                  if (sol_value>best_value)
                     {store_new_best_sol(n,gc,start,time_values_opt,
                          time_to_opt,start_time,psol,pgr);
                      best_value=sol_value;
                     }
                  continue;
                 }
              acc_val=random(seed,coef);
              db=((double)dif)/t;db=exp(db);
              if (acc_val>db) continue;
              update_temp_solution(n,k,m,pinst,psol,pgr);
              sol_value+=dif;
             }
           else
             {size=1;
              while (size>0)
                 {k=random(seed,coef)*n+1;
                  size=(pgr+(psol+k)->s)->sz;
                  if (size>(pgr+(psol+k)->s)->mins) break;
                 }
              kg=(psol+k)->s;j=kg;
              while (j==kg)
                 {j=random(seed,coef)*gc+1;
                  if (j!=kg && (pgr+j)->sz>=(pgr+j)->maxs) j=kg;
                 }
              dif=con(k,j)-con(k,kg);
              if (dif>=0)
                 {update_temp_solution(n,k,-j,pinst,psol,pgr);
                  sol_value+=dif;
                  if (sol_value>best_value)
                     {store_new_best_sol(n,gc,start,time_values_opt,
                          time_to_opt,start_time,psol,pgr);
                      best_value=sol_value;
                     }
                  continue;
                 }
              acc_val=random(seed,coef);
              db=((double)dif)/t;db=exp(db);
              if (acc_val>db) continue;
              update_temp_solution(n,k,-j,pinst,psol,pgr);
              sol_value+=dif;
             }
         } // for s
      if (start>1) {
          end_time=clock();
          elapsed_time=(long)(end_time-start_time)/CLK_TCK;
          if (elapsed_time>=time_limit) {*stop_cond=1;break;}
         }
      t=coolf*t;
     } // for r
  return best_value;
 }


 double sa_multistart_mode(int n,int gc,int rep_coef,
      int pair_exch_count,long time_limit,double min_temp,
      double coolf,double rel_prob,int *time_values_opt,
      double *time_to_opt,double *seed,clock_t start_time,
      Instance *pinst,Solution *psol,Group *pgr)
 {int stop_cond=0;
  int t_val_count;
  int rep_count;
  long start=0;
  double best_value;
  double max_temp;
  double coef;

  coef=2048;coef*=1024;coef*=1024;coef-=1;
  (psol+1)->perf=0;

// Random grouping is generated
  best_value=start_sol(n,gc,coef,seed,pinst,psol,pgr);
  store_new_best_sol(n,gc,start,time_values_opt,time_to_opt,
      start_time,psol,pgr);
  (psol+1)->perf=0;
  max_temp=set_max_temperature(n,pair_exch_count,coef,seed,pinst,psol);
  con(0,0)=max_temp;
  (psol+8)->perf=t_val_count=
      (log(min_temp)-log(max_temp))/log(coolf);
  rep_count=rep_coef*n;

  while (stop_cond==0)
     {start++;
      best_value=simulated_annealing(n,gc,t_val_count,rep_count,
          best_value,max_temp,start,time_limit,coolf,rel_prob,
          coef,&stop_cond,time_values_opt,time_to_opt,seed,
          start_time,pinst,psol,pgr);
     }
// Storing the total number of simulated annealing starts executed
  (psol+3)->perf=start;
  return best_value;
 }





// 'grouping' is the SA procedure for creating maximally diverse groups from a given set of elements
 void grouping_sa(char *in_file_name,char *out_file_name,double seed,
      long time_limit,Results *pres)
 {FILE *out,*in;
  Instance *pinst;
  Solution *psol;
  Group *pgr;

  int i,j,k,p,r;
  int n,gc;
  int minsize,maxsize=-1;
  int lo;
  int time_values[5],time_values_opt[5];
  double best_value;
  double we;
  double time_in_seconds,time_to_opt;
  double rel_prob;
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

  if (minsize>=maxsize) rel_prob=0; else rel_prob=RELOCATING_PROBAB;
  start_time=clock();
// Running simulated annealing procedure
  best_value=sa_multistart_mode(n,gc,REPETITION_COEF,SAMPLE_SIZE,
      time_limit,MIN_TEMPERATURE,COOLING_FACTOR,rel_prob,
      time_values_opt,&time_to_opt,&seed,start_time,pinst,psol,pgr);
  time_in_seconds=take_time(time_values,start_time);
  check_best_solution(out,n,gc,best_value,pinst,psol,pgr);

// Printing parameters, data characteristics and some statistics 
// regarding the performance of the algorithm
  fprintf(out,"   parameters:                                 \n");
  fprintf(out,"      COOLING_FACTOR                         = %8lf\n",
      COOLING_FACTOR);
  fprintf(out,"      MIN_TEMPERATURE                        = %8lf\n",
      MIN_TEMPERATURE);
  fprintf(out,"      maximum temperature                    = %8lf\n",
      con(0,0));
  fprintf(out,"      count of temperature values            = %8ld\n",
      (psol+8)->perf);
  fprintf(out,"      REPETITION_COEF                        = %5d\n",
      REPETITION_COEF);
  fprintf(out,"      SAMPLE_SIZE                            = %5d\n",
      SAMPLE_SIZE);
  fprintf(out,"      RELOCATING_PROBAB                      = %8lf\n",
      RELOCATING_PROBAB);
  fprintf(out,"      seed for random number generator       = %8lf\n",
      seed_saved);
  fprintf(out,"   number of elements                        = %5d\n",n);
  fprintf(out,"   number of groups                          = %5d\n",gc);
  fprintf(out,"   minimum group size                        = %5d\n",minsize);
  fprintf(out,"   maximum group size                        = %5d\n",maxsize);
  fprintf(out,"   time limit                                = %5ld\n",
      time_limit);
  fprintf(out,"   number of SA runs executed                = %5ld\n",
      (psol+3)->perf);
  fprintf(out,"   number of improvements during SA runs     = %5ld\n",
      (psol+1)->perf);
  fprintf(out,"   last improvement at SA run no.            = %5ld\n",
      (psol+2)->perf);
  fprintf(out,"   value found by SA                         = %8lf\n",
      best_value);
  (psol+4)->perf=(long)time_to_opt;
  fprintf(out,"   time to the best solution: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values_opt[1],time_values_opt[2],time_values_opt[3],
      time_values_opt[4],(long)time_to_opt);
  (psol+5)->perf=(long)time_in_seconds;
  fprintf(out,"   total time of SA: %d : %d : %d.%3d  (=%4ld seconds)\n",
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
  for (i=0;i<=n;i++) if ((pinst+i)->con!=NULL) free((pinst+i)->con);
  for (i=0;i<=n;i++) if ((pinst+i)->d!=NULL) free((pinst+i)->d);
  if (pinst!=NULL) free(pinst);
// Closing input and output files
  fclose(out);fclose(in);
 }
