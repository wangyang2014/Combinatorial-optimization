/* Program: an implementation of the iterated tabu search (ITS) algorithm 
   for the maximally diverse grouping problem. 
   Author: Gintaras Palubeckis
   Date: 2011-12-16
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
       characteristics. The structure is defined in file 'grouping_its.h'. 
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
     grouping_its(in_file_name,out_file_name,seed,20,pres);
       or just
     char in_file_name[80],out_file_name[80];
     double seed=1000;
     strcpy(in_file_name,"D:\\data\\RanReal_n012_ds_06.txt");
     strcpy(out_file_name,"D:\\temp\\RanReal_n012_ds_06.res");
     grouping_its(in_file_name,out_file_name,seed,20,NULL);
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
#include "grouping_its.h"



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
          "!!! group size violation in ITS solution: %4d   %4d\n",
          k,(pgr+k)->bsz);
      exit(1);
     }
  for (k=1;k<=gc;k++)
     {count=0;
      for (i=1;i<=n;i++) if ((psol+i)->bs==k) count++;
      if (count<(pgr+k)->mins || count>(pgr+k)->maxs)
         {fprintf(out,
          "!!! group size violation in ITS solution (2): %4d   %4d\n",
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
          "!!! some discrepancy in ITS solution values (1): %10lf   %10lf\n",
          best_value,value_obj);
      exit(1);
     }
  for (k=1;k<=gc;k++)
  for (i=1;i<(pgr+k)->bsz;i++) for (j=i+1;j<=(pgr+k)->bsz;j++)
      value_groups+=(d(bsol(k,i),bsol(k,j)));
  if (value_groups<best_value-0.00001 || value_groups>best_value+0.00001)
     {fprintf(out,
          "!!! some discrepancy in ITS solution values (2): %10lf   %10lf\n",
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


 double get_perturb_solution(int n,int gc,int perturb_count,
      int cand_list_size,int adv_use_rel,double sol_value,double coef,
      double *seed,Instance *pinst,Solution *psol,Group *pgr,
      Move *pmov)
 {int i,j,j2,k,m;
  int kg,mg;
  int moved=0;
  int ind,minind,mind,minbuck,cand_count;
  int rel=0;
  int buck_count;
  int saved=0;
  double minval,mval,dif;
  double rnum;

  buck_count=cand_list_size/BUCKET_SIZE;
  if (buck_count*BUCKET_SIZE<cand_list_size)
     {buck_count++;saved=(pmov+buck_count)->bu;
      (pmov+buck_count)->bu=cand_list_size;
     }
// sel==1 for those elements which have been moved from one group to another one
  for (i=1;i<=n;i++) (psol+i)->sel=0;

  while (moved<perturb_count)
     {cand_count=0;minval=POS_LARGE;
// Constructing a set of pairs "element-element" and "element-group" 
// as best candidates for performing a move.
// These pairs are (pmov+i)->v1, (pmov+i)->v2.
// The exchange (or relocation) cost is (pmov+i)->cr. 
// Clearly, it should be maximized.
      if (adv_use_rel==1)
         {rnum=random(seed,coef);
          if (rnum<=REL_PROB) rel=1; else rel=0;
         }
      for (k=1;k<=n;k++)
         {if ((psol+k)->sel>0) continue;
          kg=(psol+k)->s;
          if (rel>0) goto lab1;
          for (m=k+1;m<=n;m++)
             {if ((psol+m)->s==kg || (psol+m)->sel>0) continue;
              mg=(psol+m)->s;
              dif=con(k,mg)-con(k,kg)+con(m,kg)-con(m,mg)-2*d(k,m);
              if (cand_count<cand_list_size)
                 {cand_count++;
                  (pmov+cand_count)->v1=k;(pmov+cand_count)->v2=m;
                  (pmov+cand_count)->cr=dif;
                  if (cand_count<cand_list_size) continue;
                  for (i=1;i<=buck_count;i++)
                     {mind=(pmov+i)->bl;mval=(pmov+mind)->cr;
                      for (j=mind+1;j<=(pmov+i)->bu;j++)
                          if ((pmov+j)->cr<mval)
                             {mval=(pmov+j)->cr;mind=j;}
                      (pmov+i)->minv=mval;(pmov+i)->ind=mind;
                      if (mval<minval) {minval=mval;minbuck=i;}
                     } // for i
                 }
               else if (dif>minval)
                 {minind=(pmov+minbuck)->ind;
                  (pmov+minind)->v1=k;(pmov+minind)->v2=m;
                  (pmov+minind)->cr=dif;
                  mind=(pmov+minbuck)->bl;mval=(pmov+mind)->cr;
                  for (j=mind+1;j<=(pmov+minbuck)->bu;j++)
                      if ((pmov+j)->cr<mval) {mval=(pmov+j)->cr;mind=j;}
                  (pmov+minbuck)->minv=mval;(pmov+minbuck)->ind=mind;
                  minval=(pmov+1)->minv;minbuck=1;
                  for (i=2;i<=buck_count;i++)
                      if ((pmov+i)->minv<minval) {minval=(pmov+i)->minv;minbuck=i;}
                 }
             } // for m
     lab1:if ((pgr+kg)->sz<=(pgr+kg)->mins) continue;
          for (j=1;j<=gc;j++)
             {if (j==kg || (pgr+j)->sz>=(pgr+j)->maxs) continue;
              dif=con(k,j)-con(k,kg);
              if (cand_count<cand_list_size)
                 {cand_count++;
                  (pmov+cand_count)->v1=k;(pmov+cand_count)->v2=-j;
                  (pmov+cand_count)->cr=dif;
                  if (cand_count<cand_list_size) continue;
                  for (i=1;i<=buck_count;i++)
                     {mind=(pmov+i)->bl;mval=(pmov+mind)->cr;
                      for (j2=mind+1;j2<=(pmov+i)->bu;j2++)
                          if ((pmov+j2)->cr<mval)
                             {mval=(pmov+j2)->cr;mind=j2;}
                      (pmov+i)->minv=mval;(pmov+i)->ind=mind;
                      if (mval<minval) {minval=mval;minbuck=i;}
                     } // for i
                 }
               else if (dif>minval)
                 {minind=(pmov+minbuck)->ind;
                  (pmov+minind)->v1=k;(pmov+minind)->v2=-j;
                  (pmov+minind)->cr=dif;
                  mind=(pmov+minbuck)->bl;mval=(pmov+mind)->cr;
                  for (j2=mind+1;j2<=(pmov+minbuck)->bu;j2++)
                      if ((pmov+j2)->cr<mval) {mval=(pmov+j2)->cr;mind=j2;}
                  (pmov+minbuck)->minv=mval;(pmov+minbuck)->ind=mind;
                  minval=(pmov+1)->minv;minbuck=1;
                  for (i=2;i<=buck_count;i++)
                      if ((pmov+i)->minv<minval) {minval=(pmov+i)->minv;minbuck=i;}
                 }
             } // for j
         } // for k
      if (cand_count<=0) {if (rel==0) break; else continue;}
// Randomly and uniformly selecting a move from the candidate list
      ind=random(seed,coef)*cand_count+1;
      k=(pmov+ind)->v1;m=(pmov+ind)->v2;
      (psol+k)->sel=1;
      if (m>0) {(psol+m)->sel=1;moved+=2;} else moved++;
// Performing exchange (or relocation)
      sol_value+=(pmov+ind)->cr;
      update_temp_solution(n,k,m,pinst,psol,pgr);
      if (n<=moved) break;
     } // while
  if (saved>0) (pmov+buck_count)->bu=saved;
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


 double tabu_search(int n,int gc,int keep_tabu_time,int start,
      double sol_value,int it_bound,long time_limit,double coef,
      double *best_value,int *stop_cond,int *time_values_opt,
      double *time_to_opt,double *seed,clock_t start_time,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i,j,k,m;
  int ind1,ind2,imp;
  int kg,mg,prevgr;
  int count;
  int tsize=0;
  int it=0;
  double dif,val_dif;
  double temp_best;
  double elapsed_time;
  double dr_numb;
  clock_t end;

// Tabu tenure values are set to 0
  for (i=1;i<=n;i++)
     {for (j=1;j<=n;j++) texch(i,j)=0;
      for (j=1;j<=gc;j++) trel(i,j)=0;
     }
  while (it<it_bound)
     {it++;imp=0;
      ind1=ind2=-1;temp_best=NEG_LARGE;val_dif=*best_value-sol_value+0.00001;
      for (k=1;k<=n;k++)
         {kg=(psol+k)->s;
// Looking for the best pairwise exchange of elements 
// (currently best elements are denoted by ind1, ind2)
          for (m=k+1;m<=n;m++)
             {if ((psol+m)->s==kg) continue;
              mg=(psol+m)->s;
              dif=con(k,mg)-con(k,kg)+con(m,kg)-con(m,mg)-2*d(k,m);
// If the condition below is true, then the new best solution is found
              if (dif>val_dif)
                 {if (imp==0)
                     {temp_best=dif;ind1=k;ind2=m;imp=1;}
                   else
                     {imp++;dr_numb=random(seed,coef);
                      if (dr_numb<=1./(double)imp)
                         {temp_best=dif;ind1=k;ind2=m;}
                     }
                 }
// Pairs, both elements of which are in the tabu list, are bypassed
              if (imp>0 || texch(k,m)>0) continue;
              if (dif>temp_best+0.00001)
                 {temp_best=dif;ind1=k;ind2=m;count=1;}
                else if (dif>=temp_best-0.00001)
                 {count++;dr_numb=random(seed,coef);
                  if (dr_numb<=1./(double)count) {ind1=k;ind2=m;}
                 }
             } // for m
          if ((pgr+kg)->sz<=(pgr+kg)->mins) continue;
// Evaluating relocations of the element (when the element is moved 
// from its current group to another one). 
// Now (ind1, -ind2) denotes a pair (element, its new group).
// The minus sign is used to make a distinction from the case (element, element)
          for (j=1;j<=gc;j++)
             {if (j==kg || (pgr+j)->sz>=(pgr+j)->maxs) continue;
              dif=con(k,j)-con(k,kg);
              if (dif>val_dif)
                 {if (imp==0)
                     {temp_best=dif;ind1=k;ind2=-j;imp=1;}
                   else
                     {imp++;dr_numb=random(seed,coef);
                      if (dr_numb<=1./(double)imp)
                         {temp_best=dif;ind1=k;ind2=-j;}
                     }
                 }
              if (imp>0 || trel(k,j)>0) continue;
              if (dif>temp_best+0.00001)
                 {temp_best=dif;ind1=k;ind2=-j;count=1;}
                else if (dif>=temp_best-0.00001)
                 {count++;dr_numb=random(seed,coef);
                  if (dr_numb<=1./(double)count) {ind1=k;ind2=-j;}
                 }
             } // for j
         } // for k
// Assigning selected elements (or one element) to new groups (or group)
      prevgr=(psol+ind1)->s;
      update_temp_solution(n,ind1,ind2,pinst,psol,pgr);
// Solution value is updated
      sol_value+=temp_best;
// Once a new best solution is found, this solution is submitted to the local search procedure 
// for possible futher improvement
      if (imp>0)
         {sol_value=local_search(n,gc,sol_value,coef,seed,pinst,psol,pgr);
          store_new_best_sol(n,gc,start,time_values_opt,time_to_opt,start_time,psol,pgr);
          *best_value=sol_value;
         }
// Updating tabu tenure values.
      if (tsize>=keep_tabu_time)
         {k=(psol+1)->t1;m=(psol+1)->t2;
          if (m>0) texch(k,m)=0; else trel(k,-m)=0;
          for (i=1;i<tsize;i++)
             {(psol+i)->t1=(psol+i+1)->t1;(psol+i)->t2=(psol+i+1)->t2;}
         }
       else tsize++;
      (psol+tsize)->t1=ind1;(psol+tsize)->t2=ind2;
      if (ind2>0) texch(ind1,ind2)=1; else trel(ind1,prevgr)=1;
// Checking tabu search termination rule
      end=clock();
      elapsed_time=(end-start_time)/CLK_TCK;
      if (elapsed_time>=time_limit) {*stop_cond=1;break;}
     } // while
  return sol_value;
 }


 double iterated_tabu_search(int n,int gc,int keep_tabu_time,
      int perturb_count1,int perturb_count2,int min_perturb_count,
      int cand_list_size1,int cand_list_size2,int adv_use_rel,
      int it_bound,long time_limit,int *time_values_opt,
      double *time_to_opt,double *seed1,clock_t start,
      Instance *pinst,Solution *psol,Group *pgr,Move *pmov)
 {int i,s;
  int st=1;
  int stop_cond=0;
  int perturb_count,cand_list_size;
  int loind=1,hiind,count=0;
  double best_value,sol_value;
  double sum,max_val;
  double seed2,seed3,coef;

  coef=2048;coef*=1024;coef*=1024;coef-=1;
  seed2=2*(*seed1);seed3=3*(*seed1);
  if (adv_use_rel>0)
     {for (i=2,max_val=-1;i<=n;i++) if (d(1,i)>max_val) {max_val=d(1,i);s=i;}
      for (i=2;i<=n;i++)
         {if (i==s) continue;
          sum=d(1,i)+d(s,i);
          if (sum<max_val-0.00001) {adv_use_rel=0;break;}
         }
     }
  while (loind<=cand_list_size2)
     {hiind=loind+BUCKET_SIZE-1;
      if (hiind>cand_list_size2) hiind=cand_list_size2;
      count++;
      (pmov+count)->bl=loind;(pmov+count)->bu=hiind;
      loind=hiind+1;
     }

// Generating initial solution randomly
  best_value=sol_value=random_start(n,gc,coef,seed1,pinst,psol,pgr);
  store_new_best_sol(n,gc,0,time_values_opt,time_to_opt,start,psol,pgr);
  (psol+1)->perf=0;
  while (stop_cond==0)
// Invocation of tabu search
     {sol_value=tabu_search(n,gc,keep_tabu_time,st,sol_value,it_bound,
          time_limit,coef,&best_value,&stop_cond,time_values_opt,
          time_to_opt,&seed2,start,pinst,psol,pgr);
      if (stop_cond>0) break;
      st++;
      perturb_count=
          random(&seed2,coef)*(perturb_count2-perturb_count1+1);
      perturb_count+=perturb_count1;
// i stands for the number of elements for which assignement to groups will be changed
      if (perturb_count<=min_perturb_count) i=perturb_count;
       else 
         {i=random(&seed2,coef)*(perturb_count-min_perturb_count+1);
          i+=min_perturb_count;
         }
      cand_list_size=
          random(&seed2,coef)*(cand_list_size2-cand_list_size1+1);
      cand_list_size+=cand_list_size1;
// Perturbation of the current solution
      sol_value=get_perturb_solution(n,gc,i,cand_list_size,adv_use_rel,
          sol_value,coef,&seed3,pinst,psol,pgr,pmov);
     }
// Storing the total number of tabu search starts executed
  (psol+3)->perf=st;
  return best_value;
 }





// 'grouping' is the ITS procedure for creating maximally diverse groups 
// from a given set of elements
 void grouping_its(char *in_file_name,char *out_file_name,
      double seed,long time_limit,Results *pres)
 {FILE *out,*in;
  Instance *pinst;
  Solution *psol;
  Group *pgr;
  Move *pmov;

  int i,j,k,p,r;
  int n,gc;
  int minsize,maxsize=-1;
  int lo;
  int keep_tabu_time;
  int perturb_count1,perturb_count2;
  int it_bound;
  int adv_use_rel=0;
  int time_values[5],time_values_opt[5];
  double best_value;
  double we;
  double time_in_seconds,time_to_opt;
  double seed_saved;
  clock_t start;
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
  for (i=0;i<=n;i++) ALI((psol+i)->texch,n+1)
  for (i=0;i<=n;i++) ALI((psol+i)->trel,gc+1)
  ALS(pmov,Move,LIST_SIZE2+1)
// Allocation of core memory for an array of Group structure objects 
// that contains data used in the construction of groups
  ALS(pgr,Group,gc+1)
  minsize=n+1;
  for (i=1;i<=gc;i++)
     {fscanf(in,"%d %d",&p,&r);
      (pgr+i)->mins=p;(pgr+i)->maxs=r;
      if (p<r) adv_use_rel=1;
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

// Preparation of ITS parameters
  keep_tabu_time=TABU_TIME;i=n/TABU_COEF;
  if (i<keep_tabu_time) keep_tabu_time=i;
  perturb_count1=n*PER_COEF1;perturb_count2=n*PER_COEF2;
  if (n>=PROBL_SIZE_CONST) it_bound=ITER_COUNT_BIGGER;
   else it_bound=ITER_COUNT_SMALLER;

  start=clock();
// Running iterated tabu search
  best_value=iterated_tabu_search(n,gc,keep_tabu_time,
      perturb_count1,perturb_count2,MIN_PER_COUNT,LIST_SIZE1,
      LIST_SIZE2,adv_use_rel,it_bound,time_limit,time_values_opt,
      &time_to_opt,&seed,start,pinst,psol,pgr,pmov);
  time_in_seconds=take_time(time_values,start);
  check_best_solution(out,n,gc,best_value,pinst,psol,pgr);

// Printing parameters, data characteristics and some statistics 
// regarding the performance of the algorithm
  fprintf(out,"   parameters:                                 \n");
  fprintf(out,"      TABU_TIME                              = %5d\n",
      TABU_TIME);
  fprintf(out,"      TABU_COEF                              = %5d\n",
      TABU_COEF);
  fprintf(out,"      PER_COEF1                              = %6lf\n",
      PER_COEF1);
  fprintf(out,"      PER_COEF2                              = %6lf\n",
      PER_COEF2);
  fprintf(out,"      MIN_PER_COUNT                          = %5d\n",
      MIN_PER_COUNT);
  fprintf(out,"      LIST_SIZE1                             = %5d\n",
      LIST_SIZE1);
  fprintf(out,"      LIST_SIZE2                             = %5d\n",
      LIST_SIZE2);
  fprintf(out,"      REL_PROB                               = %5lf\n",
      REL_PROB);
  fprintf(out,"      number of iterations per TS start      = %5d\n",
      it_bound);
  fprintf(out,"      seed for random number generator       = %8lf\n",
      seed_saved);
  fprintf(out,"   number of elements                        = %5d\n",n);
  fprintf(out,"   number of groups                          = %5d\n",gc);
  fprintf(out,"   minimum group size                        = %5d\n",minsize);
  fprintf(out,"   maximum group size                        = %5d\n",maxsize);
  fprintf(out,"   time limit                                = %5ld\n",
      time_limit);
  fprintf(out,"   number of TS starts executed              = %4ld\n",
      (psol+3)->perf);
  fprintf(out,"   number of improvements during ITS         = %4ld\n",
      (psol+1)->perf);
  fprintf(out,"   last improvement at TS start no.          = %3ld\n",
      (psol+2)->perf);
  fprintf(out,"   value found by ITS                        = %8lf\n",
      best_value);
  (psol+4)->perf=(long)time_to_opt;
  fprintf(out,"   time to the best solution: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values_opt[1],time_values_opt[2],time_values_opt[3],
      time_values_opt[4],(long)time_to_opt);
  (psol+5)->perf=(long)time_in_seconds;
  fprintf(out,"   total time of ITS: %d : %d : %d.%3d  (=%4ld seconds)\n",
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
  if (pmov!=NULL) free(pmov);
  for (i=0;i<=n;i++) if ((psol+i)->trel!=NULL) free((psol+i)->trel);
  for (i=0;i<=n;i++) if ((psol+i)->texch!=NULL) free((psol+i)->texch);
  if (psol!=NULL) free(psol);
  for (i=0;i<=n;i++) if ((pinst+i)->con!=NULL) free((pinst+i)->con);
  for (i=0;i<=n;i++) if ((pinst+i)->d!=NULL) free((pinst+i)->d);
  if (pinst!=NULL) free(pinst);
// Closing input and output files
  fclose(out);fclose(in);
 }
