/* Program: an implementation of the genetic algorithm 
   for the maximally diverse grouping problem. 
   Author: Gintaras Palubeckis
   Date: 2011-07-08
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
       characteristics. The structure is defined in file 'grouping_ga.h'. 
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
     grouping_ga(in_file_name,out_file_name,seed,20,pres);
       or just
     char in_file_name[80],out_file_name[80];
     double seed=1000;
     strcpy(in_file_name,"D:\\data\\RanReal_n012_ds_06.txt");
     strcpy(out_file_name,"D:\\temp\\RanReal_n012_ds_06.res");
     grouping_ga(in_file_name,out_file_name,seed,20,NULL);
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
#include "grouping_ga.h"



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


 void sort_int_pairs_decr(int *weaux,int *we,int *caux,
      int *c,long s_ind,long e_ind)
 {int i,j,k,l;
  int m_ind;
  int e;
  int gap;

// This code implements the well-known merge algorithm for the sorting problem.
// The array submitted for sorting is c(.). 
// It is sorted in nonincreasing we(.) order
  gap=e_ind-s_ind;
  if (gap<7)
     {for (i=s_ind;i<e_ind;i++)
          for (j=i;j>s_ind && *(we+j-1)<*(we+j);j--)
             {e=*(we+j-1);*(we+j-1)=*(we+j);*(we+j)=e;
              e=*(c+j-1);*(c+j-1)=*(c+j);*(c+j)=e;
             }
      return;
     }
  m_ind=(s_ind+e_ind)/2;
  sort_int_pairs_decr(we,weaux,c,caux,s_ind,m_ind);
  sort_int_pairs_decr(we,weaux,c,caux,m_ind,e_ind);
  if (*(weaux+m_ind-1)>=*(weaux+m_ind))
     {for (i=s_ind;i<e_ind;i++)
         {*(we+i)=*(weaux+i);*(c+i)=*(caux+i);}
      return;
     }
  for (i=s_ind,k=s_ind,l=m_ind;i<e_ind;i++)
     {if (l>=e_ind || k<m_ind && *(weaux+k)>=*(weaux+l))
         {*(we+i)=*(weaux+k);*(c+i)=*(caux+k);
          k++;
         }
       else
         {*(we+i)=*(weaux+l);*(c+i)=*(caux+l);
          l++;
         }
     }
 }


 void check_best_solution(FILE *out,int n,int gc,double best_value,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i,j,k;
  int count;
  double value_obj=0,value_groups=0;

  for (k=1;k<=gc;k++) if ((pgr+k)->bsz<(pgr+k)->mins || (pgr+k)->bsz>(pgr+k)->maxs)
     {fprintf(out,
          "!!! group size violation in GA solution: %4d   %4d\n",
          k,(pgr+k)->bsz);
      exit(1);
     }
  for (k=1;k<=gc;k++)
     {count=0;
      for (i=1;i<=n;i++) if ((psol+i)->bs==k) count++;
      if (count<(pgr+k)->mins || count>(pgr+k)->maxs)
         {fprintf(out,
          "!!! group size violation in GA solution (2): %4d   %4d\n",
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
          "!!! some discrepancy in GA solution values (1): %10lf   %10lf\n",
          best_value,value_obj);
      exit(1);
     }
  for (k=1;k<=gc;k++)
  for (i=1;i<(pgr+k)->bsz;i++) for (j=i+1;j<=(pgr+k)->bsz;j++)
      value_groups+=(d(bsol(k,i),bsol(k,j)));
  if (value_groups<best_value-0.00001 || value_groups>best_value+0.00001)
     {fprintf(out,
          "!!! some discrepancy in GA solution values (2): %10lf   %10lf\n",
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


void preparation_for_local_search(int n,int gc,Instance *pinst,
      Solution *psol)
 {int i,j,k,m;
  int kg,mg;

  for (i=1;i<=n;i++) for (j=1;j<=gc;j++) con(i,j)=0;
  for (k=1;k<n;k++)
     {kg=(psol+k)->s;
      for (m=k+1;m<=n;m++)
         {mg=(psol+m)->s;
          (con(k,mg))+=(d(k,m));(con(m,kg))+=(d(k,m));
         }
     }
 }


 double local_search(int n,int gc,int prepar_req,double best_value,
      double coef,double *seed,Instance *pinst,Solution *psol,
      Group *pgr)
 {int i,j,j0,k,k0,m,m0;
  int kg,mg;
  int ind1,ind2;
  double dif;
  double temp_best=1;
  double sol_value;

  sol_value=best_value;
  if (prepar_req>0) preparation_for_local_search(n,gc,pinst,psol);
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


 void store_new_best_sol(int n,int gc,double new_value,
      long offspr_no,int *time_values_opt,
      double *best_value,double *time_to_opt,
      clock_t start_time,Solution *psol,Group *pgr)
 {int i,j;

  *best_value=new_value;
  for (j=1;j<=gc;j++) {(pgr+j)->bsz=(pgr+j)->sz;(pgr+j)->c=0;}
  for (i=1;i<=n;i++)
     {(psol+i)->bs=(psol+i)->s;
      j=(psol+i)->s;((pgr+j)->c)++;bsol(j,(pgr+j)->c)=i;
     }
  if (offspr_no>0) {((psol+1)->perf)++;(psol+2)->perf=offspr_no;}
  *time_to_opt=take_time(time_values_opt,start_time);
 }


 void append_to_population(int n,int gc,int memb_count,double sol_value,
      int *time_values_opt,double *best_value,double *time_to_opt,
      clock_t start_time,Solution *psol,Group *pgr)
 {int i,j;
  int ln;

  if (memb_count==0 || sol_value>*best_value)
      store_new_best_sol(n,gc,sol_value,0,time_values_opt,
          best_value,time_to_opt,start_time,psol,pgr);
  ln=memb_count+1;
  for (i=1;i<=n;i++) pop(ln,i)=(psol+i)->s;
  for (j=1;j<=gc;j++) popsz(j,ln)=(pgr+j)->sz;
  pfv(ln)=sol_value;
  if (ln==1 || sol_value<pfv(0)) {pfv(0)=sol_value;psol->r=ln;}
 }


 void sort_by_max_size(int gc,Group *pgr)
 {int i;

  for (i=1;i<=gc;i++)
     {oob(0,i)=i;oob(1,i)=i;
      p(0,i)=(pgr+i)->maxs;p(1,i)=p(0,i);
     }
  sort_int_pairs_decr((pgr+1)->p,pgr->p,(pgr+1)->oob,pgr->oob,1,gc+1);
 }


 int init_population(int n,int gc,int psize_req,double coef,
      int *time_values_opt,double *best_value,double *time_to_opt,
      double *seed,clock_t start_time,Instance *pinst,
      Solution *psol,Group *pgr)
 {int memb_count=0;
  double memb_value;

  while (memb_count<psize_req)
     {memb_value=random_start(n,gc,coef,seed,pinst,psol,pgr);
      memb_value=local_search(n,gc,0,memb_value,coef,seed,pinst,psol,pgr);
      append_to_population(n,gc,memb_count,memb_value,
          time_values_opt,best_value,time_to_opt,start_time,psol,pgr);
      memb_count++;
     }
  return memb_count;
 }


 void select_parents(int psize,double coef,double *seed,
      int *pr1,int *pr2)
 {int i;

  *pr1=random(seed,coef)*psize+1;
  i=random(seed,coef)*(psize-1)+1;
  if (i>=*pr1) i++;
  *pr2=i;
 }


 int remove_excessive(int r,int tot_count,double coef,
      double *seed,Solution *psol,Group *pgr)
 {int i,k,s;
  int size;

  size=(pgr+r)->oc;
  for (i=1;i<=size;i++)
     {k=random(seed,coef)*(size-i+1);k+=i;
      s=oob(r,i);oob(r,i)=oob(r,k);oob(r,k)=s;
     }
  for (i=(pgr+r)->maxs+1;i<=size;i++)
     {(psol+oob(r,i))->sel=0;tot_count--;}
  return tot_count;
 }


 int make_free(int gc,int def,int tot_count,double coef,
      double *seed,Solution *psol,Group *pgr)
 {int i,j,k,m,s;
  int count=0;
  int grp;
  int size;

  for (i=1;i<=gc;i++)
     {m=(pgr+i)->oc-(pgr+i)->mins;
      if (m<=0) continue;
      for (j=1;j<=m;j++) {count++;(psol+count)->r=i;}
     }
  for (i=1;i<=count;i++)
     {k=random(seed,coef)*(count-i+1);k+=i;
      s=(psol+i)->r;(psol+i)->r=(psol+k)->r;(psol+k)->r=s;
     }
  for (i=1;i<=def;i++)
     {grp=(psol+i)->r;size=(pgr+grp)->oc;
      k=random(seed,coef)*size+1;
      (psol+oob(grp,k))->sel=0;tot_count--;
      oob(grp,k)=oob(grp,size);
      ((pgr+grp)->oc)--;
     }
  return tot_count;
 }


 double get_offspring(int n,int gc,int maxsize,
      int pr1,int pr2,double coef,double *seed,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i,j,k,l,r,s,t;
  int ig,grp1,grp2;
  int count;
  int sz,sz_some=-1;
  int tot_count=0,gcount=0,rem,cand=0;
  int need,take,deficit;
  double offspr_value=0;

  for (i=0;i<=maxsize;i++) (psol+i)->c=0;
  for (i=1;i<=gc;i++) (pgr+i)->c1=(pgr+i)->c2=(pgr+i)->oc=0;
  for (i=1;i<=n;i++)
     {ig=pop(pr1,i);((pgr+ig)->c1)++;
      ob1(ig,(pgr+ig)->c1)=i;
      ig=pop(pr2,i);((pgr+ig)->c2)++;
      ob2(ig,(pgr+ig)->c2)=i;
      (psol+i)->sel=0;
     }
  for (i=1;i<=gc;i++)
     {for (k=1;k<=(pgr+i)->c1;k++) (psol+ob1(i,k))->m=1;
      for (j=1;j<=gc;j++)
         {for (l=1;l<=(pgr+j)->c2;l++) ((psol+ob2(j,l))->m)++;
          for (k=1,count=0;k<=(pgr+i)->c1;k++)
              if ((psol+ob1(i,k))->m>1) count++;
          is(i,j)=count;((psol+count)->c)++;
          for (l=1;l<=(pgr+j)->c2;l++) ((psol+ob2(j,l))->m)--;
         }
      for (k=1;k<=(pgr+i)->c1;k++) (psol+ob1(i,k))->m=0;
     }

  for (i=maxsize,need=gc;i>=1;i--)
     {if ((psol+i)->c<need) {need-=(psol+i)->c;continue;}
      if ((psol+i)->c==need) {sz=i;break;}
      sz=i+1;sz_some=i;break;
     }
  for (i=1;i<=gc;i++) for (j=1;j<=gc;j++)
     {if (is(i,j)>=sz)
         {gcount++;
          (pgr+gcount)->i1=i;(pgr+gcount)->i2=j;
          is(0,gcount)=gcount;bsol(0,gcount)=gcount;
          ob1(0,gcount)=is(i,j);ob2(0,gcount)=ob1(0,gcount);
          continue;
         }
      if (is(i,j)!=sz_some) continue;
      cand++;
      p(0,cand)=cand;p(1,cand)=i;p(2,cand)=j;
     }
  if (sz_some>0)
     {for (i=1;i<=cand;i++)
         {k=random(seed,coef)*(cand-i+1);k+=i;
          s=p(0,i);p(0,i)=p(0,k);p(0,k)=s;
         }
      rem=gc-gcount;
      for (i=1;i<=rem;i++)
         {gcount++;
          (pgr+gcount)->i1=p(1,p(0,i));(pgr+gcount)->i2=p(2,p(0,i));
          is(0,gcount)=gcount;bsol(0,gcount)=gcount;
          ob1(0,gcount)=sz_some;ob2(0,gcount)=ob1(0,gcount);
         }
     }

  sort_int_pairs_decr(pgr->ob2,pgr->ob1,pgr->bsol,pgr->is,1,gc+1);
  for (r=1;r<=gc;r++)
     {p(1,r)=(pgr+r)->i1;p(2,r)=(pgr+r)->i2;}
  for (r=1;r<=gc;r++)
     {k=is(0,r);j=oob(0,r);
      (pgr+j)->i1=p(1,k);(pgr+j)->i2=p(2,k);
     }

  for (r=1;r<=gc;r++)
     {i=(pgr+r)->i1;j=(pgr+r)->i2;
      for (k=1;k<=(pgr+i)->c1;k++) (psol+ob1(i,k))->m=1;
      for (l=1;l<=(pgr+j)->c2;l++) ((psol+ob2(j,l))->m)++;
      for (k=1;k<=(pgr+i)->c1;k++)
          if ((psol+ob1(i,k))->m>1)
             {((pgr+r)->oc)++;oob(r,(pgr+r)->oc)=ob1(i,k);
              (psol+ob1(i,k))->sel=1;tot_count++;
             }
      for (l=1;l<=(pgr+j)->c2;l++) ((psol+ob2(j,l))->m)--;
      for (k=1;k<=(pgr+i)->c1;k++) (psol+ob1(i,k))->m=0;
     }

  for (r=1,need=0;r<=gc;r++)
     {if ((pgr+r)->oc<=(pgr+r)->maxs)
         {if ((pgr+r)->oc<(pgr+r)->mins) need+=((pgr+r)->mins-(pgr+r)->oc);
          continue;
         }
      tot_count=remove_excessive(r,tot_count,coef,seed,psol,pgr);
      (pgr+r)->oc=(pgr+r)->maxs;
     }
  deficit=need-n+tot_count;
  if (deficit>0)
      tot_count=make_free(gc,deficit,tot_count,coef,seed,psol,pgr);

  for (r=1;r<=gc;r++)
     {if ((pgr+r)->oc>=(pgr+r)->mins) continue;
      i=(pgr+r)->i1;j=(pgr+r)->i2;rem=0;
      for (k=1;k<=(pgr+i)->c1;k++) if ((psol+ob1(i,k))->sel==0)
         {rem++;(psol+rem)->r=ob1(i,k);}
      for (l=1;l<=(pgr+j)->c2;l++) if ((psol+ob2(j,l))->sel==0)
         {rem++;(psol+rem)->r=ob2(j,l);}
      if (rem==0) continue;
      need=(pgr+r)->mins-(pgr+r)->oc;
      if (rem<=need) take=rem;
       else
         {take=need;
          for (i=1;i<=rem;i++)
             {k=random(seed,coef)*(rem-i+1);k+=i;
              s=(psol+i)->r;(psol+i)->r=(psol+k)->r;(psol+k)->r=s;
             }
         }
      for (i=1;i<=take;i++)
         {((pgr+r)->oc)++;oob(r,(pgr+r)->oc)=(psol+i)->r;
          (psol+(psol+i)->r)->sel=1;tot_count++;
         }
     }
  if (tot_count>=n) goto lab1;

  for (i=1,rem=0;i<=n;i++)
      if ((psol+i)->sel==0) {rem++;(psol+rem)->r=i;}
  for (r=1;r<=gc;r++)
     {if ((pgr+r)->oc>=(pgr+r)->mins) continue;
      need=(pgr+r)->mins-(pgr+r)->oc;
      for (i=1;i<=need;i++)
         {k=random(seed,coef)*rem+1;
          s=(psol+k)->r;
          ((pgr+r)->oc)++;oob(r,(pgr+r)->oc)=s;
          (psol+s)->sel=1;tot_count++;
          (psol+k)->r=(psol+rem)->r;
          rem--;
         }
     }
  if (rem<=0) goto lab1;

  for (i=1;i<=rem;i++)
     {s=(psol+i)->r;grp1=pop(pr1,s);grp2=pop(pr2,s);
      for (r=1,count=0;r<=gc;r++)
         {if (grp1!=(pgr+r)->i1 && grp2!=(pgr+r)->i2) continue;
          if ((pgr+r)->oc>=(pgr+r)->maxs) continue;
          count++;(pgr+count)->cand=r;
         }
      if (count==0) continue;
      k=random(seed,coef)*count+1;
      r=(pgr+k)->cand;
      ((pgr+r)->oc)++;oob(r,(pgr+r)->oc)=s;
      (psol+s)->sel=1;tot_count++;
     }
  if (tot_count>=n) goto lab1;

  for (i=1;i<=rem;i++)
     {s=(psol+i)->r;if ((psol+s)->sel>0) continue;
      for (r=1,count=0;r<=gc;r++)
         {if ((pgr+r)->oc>=(pgr+r)->maxs) continue;
          count++;(pgr+count)->cand=r;
         }
      k=random(seed,coef)*count+1;
      r=(pgr+k)->cand;
      ((pgr+r)->oc)++;oob(r,(pgr+r)->oc)=s;
      (psol+s)->sel=1;//tot_count++;
     }

 lab1:for (r=1;r<=gc;r++)
     {(pgr+r)->sz=(pgr+r)->oc;
      for (k=1;k<=(pgr+r)->oc;k++)
         {j=oob(r,k);(psol+j)->s=r;
          for (l=k+1;l<=(pgr+r)->oc;l++)
             {t=oob(r,l);offspr_value+=(d(j,t));}
         }
     }
  return offspr_value;
 }


 int evaluate_offspring(int n,int gc,int psize,int psize_req,
      double offspr_value,long offspr_no,int *time_values_opt,
      double *best_value,double *time_to_opt,
      clock_t start_time,Solution *psol,Group *pgr)
 {int i,j;
  int ind,replaced=0;
  double minval;

  if (offspr_value>*best_value)
     {store_new_best_sol(n,gc,offspr_value,offspr_no,time_values_opt,
          best_value,time_to_opt,start_time,psol,pgr);
     }
  if (offspr_value<pfv(0) && psize>=psize_req) return psize;
  if (psize<psize_req) {psize++;ind=psize;}
   else {ind=psol->r;replaced=1;}
  for (i=1;i<=n;i++) pop(ind,i)=(psol+i)->s;
  for (j=1;j<=gc;j++) popsz(j,ind)=(pgr+j)->sz;
  pfv(ind)=offspr_value;
  if (replaced>0)
     {minval=pfv(1);ind=1;
      for (i=2;i<=psize;i++) if (pfv(i)<minval)
         {minval=pfv(i);ind=i;}
      pfv(0)=minval;psol->r=ind;
     }
   else if (offspr_value<pfv(0))
     {pfv(0)=offspr_value;psol->r=psize;}
  return psize;
 }


 double genetic_algorithm(FILE *out,int n,int gc,int maxsize,
      int psize_req,long time_limit,int *time_values_opt,
      double *time_to_opt,double *seed1,clock_t start_time,
      Instance *pinst,Solution *psol,Group *pgr)
 {int i;
  int psize;
  int pr1,pr2;
  int stop_cond=0;
  long offspr_no=0;
  double best_value,offspr_value;
  double seed2,seed3,coef;
  double elapsed_time;
  clock_t end;

  coef=2048;coef*=1024;coef*=1024;coef-=1;
  seed2=2*(*seed1);seed3=3*(*seed1);
  (psol+1)->perf=0;(psol+2)->perf=-1;
  for (i=1;i<=n;i++) (psol+i)->m=0;

  sort_by_max_size(gc,pgr);
// Initial population is generated
  psize=init_population(n,gc,psize_req,coef,time_values_opt,
      &best_value,time_to_opt,seed1,start_time,pinst,psol,pgr);
  if (psize<2)
     {fprintf(out,"!!!!!  Initial population is too small to proceed: psize = %3d\n",psize);
      exit(1);
     }
  while (stop_cond==0)
     {offspr_no++;
      select_parents(psize,coef,&seed2,&pr1,&pr2);
      offspr_value=get_offspring(n,gc,maxsize,pr1,pr2,coef,&seed3,pinst,psol,pgr);
      offspr_value=local_search(n,gc,1,offspr_value,coef,&seed3,pinst,psol,pgr);
      psize=evaluate_offspring(n,gc,psize,psize_req,offspr_value,offspr_no,
          time_values_opt,&best_value,time_to_opt,start_time,psol,pgr);
      end=clock();
      elapsed_time=(end-start_time)/CLK_TCK;
      if (elapsed_time>=time_limit) stop_cond=1;
     }
  (psol+3)->perf=offspr_no;
  return best_value;
 }





// 'grouping' is the GA procedure for creating maximally diverse groups from 
// a given set of elements
 void grouping_ga(char *in_file_name,char *out_file_name,double seed,
      long time_limit,Results *pres)
 {FILE *out,*in;
  Instance *pinst;
  Solution *psol;
  Group *pgr;

  int i,j,k,p,r;
  int n,gc;
  int minsize,maxsize=-1;
  int gr_pairs;
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
  gr_pairs=gc*gc;
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
  for (i=0;i<=n;i++) ALI((psol+i)->pop,POPULATION_SIZE+1)
  ALD(psol->pfv,POPULATION_SIZE+1)
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
  for (i=0;i<=gc;i++) ALI((pgr+i)->oob,maxsize+1)
  for (i=0;i<=gc;i++) ALI((pgr+i)->ob1,maxsize+1)
  for (i=0;i<=gc;i++) ALI((pgr+i)->ob2,maxsize+1)
  for (i=0;i<=gc;i++) ALI((pgr+i)->is,gc+1)
  for (i=0;i<=gc;i++) ALI((pgr+i)->popsz,POPULATION_SIZE+1)
  for (i=0;i<=2;i++) ALI((pgr+i)->p,gr_pairs+1)

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

  start_time=clock();
// Running genetic algorithm
  best_value=genetic_algorithm(out,n,gc,maxsize,POPULATION_SIZE,
      time_limit,time_values_opt,&time_to_opt,&seed,start_time,
      pinst,psol,pgr);
  time_in_seconds=take_time(time_values,start_time);
  check_best_solution(out,n,gc,best_value,pinst,psol,pgr);

// Printing parameters, data characteristics and some statistics 
// regarding the performance of the algorithm
  fprintf(out,"   parameters:                                 \n");
  fprintf(out,"      POPULATION_SIZE                        = %5d\n",
      POPULATION_SIZE);
  fprintf(out,"      seed for random number generator       = %10lf\n",
      seed_saved);
  fprintf(out,"   number of elements                        = %5d\n",n);
  fprintf(out,"   number of groups                          = %5d\n",gc);
  fprintf(out,"   minimum group size                        = %5d\n",minsize);
  fprintf(out,"   maximum group size                        = %5d\n",maxsize);
  fprintf(out,"   time limit                                = %8ld\n",
      time_limit);
  fprintf(out,"   number of offsprings produced             = %8ld\n",
      (psol+3)->perf);
  fprintf(out,"   number of improvements during GA          = %3ld\n",
      (psol+1)->perf);
  if ((psol+2)->perf<0)
      fprintf(out,"   best solution has been found in initial population \n");
   else
      fprintf(out,"   last improvement is due offspring no.     = %8ld\n",
          (psol+2)->perf);
  fprintf(out,"   value of the best solution                = %8lf\n",
      best_value);
  (psol+4)->perf=(long)time_to_opt;
  fprintf(out,"   time to the best solution: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values_opt[1],time_values_opt[2],time_values_opt[3],
      time_values_opt[4],(long)time_to_opt);
  (psol+5)->perf=(long)time_in_seconds;
  fprintf(out,"   total solution time: %d : %d : %d.%3d  (=%4ld seconds)\n",
      time_values[1],time_values[2],time_values[3],
      time_values[4],(long)time_in_seconds);
  fflush(out);

  if (pres!=NULL)
     {pres->value=best_value;
      pres->total_time=time_in_seconds;
      pres->time_to_best=time_to_opt;
     }

// Releasing the memory
  for (i=0;i<=2;i++) if ((pgr+i)->p!=NULL) free((pgr+i)->p);
  for (i=0;i<=gc;i++) if ((pgr+i)->popsz!=NULL) free((pgr+i)->popsz);
  for (i=0;i<=gc;i++) if ((pgr+i)->is!=NULL) free((pgr+i)->is);
  for (i=0;i<=gc;i++) if ((pgr+i)->ob2!=NULL) free((pgr+i)->ob2);
  for (i=0;i<=gc;i++) if ((pgr+i)->ob1!=NULL) free((pgr+i)->ob1);
  for (i=0;i<=gc;i++) if ((pgr+i)->oob!=NULL) free((pgr+i)->oob);
  for (i=0;i<=gc;i++) if ((pgr+i)->bsol!=NULL) free((pgr+i)->bsol);
  if (pgr!=NULL) free(pgr);
  if (psol->pfv!=NULL) free(psol->pfv);
  for (i=0;i<=n;i++) if ((psol+i)->pop!=NULL) free((psol+i)->pop);
  if (psol!=NULL) free(psol);
  for (i=0;i<=n;i++) if ((pinst+i)->con!=NULL) free((pinst+i)->con);
  for (i=0;i<=n;i++) if ((pinst+i)->d!=NULL) free((pinst+i)->d);
  if (pinst!=NULL) free(pinst);
// Closing input and output files
  fclose(out);fclose(in);
 }
