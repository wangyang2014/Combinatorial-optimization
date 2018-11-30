#include <alloc.h>
#include <process.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct
     {double value;         /* solution value                      */
      double time_to_best;  /* time to the VNS best solution, secs */
      double total_time;    /* total time, secs                    */
     }Results;

void grouping_vns(char *,char *,double,long,Results *);


 void main(int argc,char **argv)
 {FILE *out;
  Results *pres;
  char in_file_name[80];
  char out_file_name[80];
  char out_file_name_init[80];
  char summary_file_name[80];
  char bkv_string[80];
  int i,j,k;
  int count;
  long time_limit;
  double max_value=-1,best_known_value=0;
  double dif_best;
  double seeds [11]={0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 
      8000, 9000, 10000};
  char numbs[31][2]={ {'0'}, {'1'}, {'2'}, {'3'}, {'4'}, {'5'},
      {'6'}, {'7'}, {'8'}, {'9'}, {'1','0'} };
  char res_ext[4]={'.', 'r', 'e', 's'};
  double av_value=0.,av_time=0.,av_time_tot=0.;
  double dif_aver;

  if (argc<=3) {printf("  specify data, output and summary files");exit(1);}
  strcpy(in_file_name,argv[1]);
  strcpy(out_file_name_init,argv[2]);
  strcpy(summary_file_name,argv[3]);
  if (argc==5)
     {strcpy(bkv_string,argv[4]);best_known_value=atol(bkv_string);}
  pres=(Results *)calloc(1,sizeof(Results));
  if ((out=fopen(summary_file_name,"w"))==NULL)
     {printf("  fopen failed for output  %s",summary_file_name);exit(1);}

  time_limit=600;count=10;
  for (i=1;i<=count;i++)
     {strcpy(out_file_name,out_file_name_init);
      k=strlen(out_file_name);
      out_file_name[k]='_';k++;
      if (i<10) out_file_name[k]=numbs[i][0];
       else {out_file_name[k]=numbs[i][0];k++;out_file_name[k]=numbs[i][1];}
      for (j=0;j<=3;j++) out_file_name[k+j+1]=res_ext[j];
      out_file_name[k+j+1]='\0';
      grouping_vns(in_file_name,out_file_name,seeds[i],time_limit,pres);
      fprintf(out," %8lf       %8lf       %8lf\n",pres->value,
          pres->time_to_best,pres->total_time);
      av_value+=pres->value;
      av_time+=pres->time_to_best;av_time_tot+=pres->total_time;
      if (pres->value>max_value) max_value=pres->value;
     }
  av_value/=count;av_time/=count;av_time_tot/=count;
  fprintf(out,"----------------------------------------\n");
  fprintf(out,"%11.3lf    %11.3lf    %11.3lf\n",av_value,av_time,av_time_tot);
  fprintf(out,"maximum value over all runs = %8lf \n",max_value);
  if (argc==5)
     {dif_best=best_known_value-max_value;dif_aver=best_known_value-av_value;
      fprintf(out,"best known value = %8lf \n",best_known_value);
      fprintf(out,"dif_best = %8lf  dif_aver = %8lf\n",dif_best,dif_aver);
     }
  free(pres);
  fclose(out);
 }
