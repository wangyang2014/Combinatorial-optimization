/******************************************************************************************/
/***************************    0. Head Files needed     **********************************/
/******************************************************************************************/
#define WIN32_LEAN_AND_MEAN		
#include <cstdlib>
#include <cstdlib>
#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
#include <ctype.h>
#include <vector>
using namespace std;
/******************************************************************************************/
/********************    1. Data Structure and Global Variables    ************************/
/******************************************************************************************/
typedef struct Solution{
	int *p ;
	double *SizeG; 
	double cost ;
}Solution; 

typedef struct Neighborhood{
	int  type ;
	int  v ;
    int  g ;
    int  x ;
    int  y ;	
}Neighborhood;

typedef struct Pair{
    int i;
    int g; 
}Pair;

char * File_Name;
char * Output_File_Name;
char * Solution_File;

int N, K;  // node number and group number
double f, f_best; 
Solution CS, NS, GS, OS;  
Neighborhood *Neighbors;
Pair *PairSet;
double total_time, starting_time, Time_limit;
int * p;      // patition array for each vertex
int * bestp ; // patition array for each vertex 
double * SizeG ; 
int **TabuTenure; 
double ** Delta_Matrix; // incremental matrix 
double ** Delta_Matrix1; 
double ** Delta; 
double ** Delta1; 
double ** D; // distance matrix between elements
double ** DT;
double * LB; // Lower bound of the number of elements in the group i 
double * UB; // Upper bound of the number of elements in the group i 
double * w;  // node weights
/******************************************************************************************/
/********************    2. Inputing Data and Assign Memeries   ***************************/
/******************************************************************************************/
void inputing()
{
    int i,j,k; 
    int x1, x2;
    float d;
	ifstream FIC;
	FIC.open(File_Name);
	if ( FIC.fail() )
	{
		cout << "### Erreur open, File_Name " << File_Name << endl;
		exit(0);
	}
	FIC >> N >> K;
	char StrReading[100];
	FIC >> StrReading;
	if ( FIC.eof() )
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}	
	if ( strcmp(StrReading, "ds" )==0 || strcmp(StrReading, "ss" )==0 )
	{
         LB = new double [K];
         UB = new double [K]; 
         for(i=0;i<K;i++) { FIC >> LB[i];  FIC >> UB[i]; }
		// for(i=0;i<K;i++)printf("%lf %lf\n",LB[i],UB[i]);
    }
    FIC >> StrReading;
	if ( FIC.eof() )
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}	
	if ( strcmp(StrReading, "W" )==0)
	{
		 w  = new double [N];
		 for(i=0;i<N;i++) FIC >> w[i];
    }
//	for(i=0;i<N;i++)printf("%lf \n",w[i]);
    D = new double * [N];
    for(i=0;i<N;i++) D[i] = new double [N];
  //  for(i=0;i<N;i++)
  //    for(j=0;j<N;j++) D[i][j] = 0.0;  
      
    DT = new double * [N];
    for(i=0;i<N;i++) DT[i] = new double [N]; 
    
  //  for(i=0;i<N;i++)
  //    for(j=0;j<N;j++) DT[i][j] = 0.0; 
  
	while ( ! FIC.eof() )
	{
			FIC >> x1 >> x2 >> d;
			//cout << x1 <<"  "<< x2 <<"  "<<d<<" "<< endl;
			if ( x1<0 || x2<0 || x1>=N || x2 >=N )
			{
				cout << "### Error of node : x1="
					<< x1 << ", x2=" << x2 << endl;
				exit(0);
			}
			if(x1 != x2)
			{
				
                D[x2][x1] = d; 
                D[x1][x2] = D[x2][x1] ;
                DT[x2][x1] = 2.0*d;
                DT[x1][x2] = DT[x2][x1];
			}
	
     }
    for(i=0;i<N;i++) DT[i][i] = D[i][i] = 0.0; 
    FIC.close();
    //for(i=0;i<N;i++)
    //for(j=i+1;j<N;j++) printf("%lf\n",D[i][j]); 
 }
void inputing1()
{
    int i,j,k; 
    int x1, x2;
    float d;
    double U;
	ifstream FIC;
	FIC.open(File_Name);
	if ( FIC.fail() )
	{
		cout << "### Erreur open, File_Name " << File_Name << endl;
		exit(0);
	}
	FIC >> N >> K;
	if ( FIC.eof() )
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}	

    LB = new double [K];
    UB = new double [K]; 
    
    FIC >> U; 
    for(i=0;i<K;i++){ LB[i] = 0; UB[i] = U; }
    w  = new double [N];
    for(i=0;i<N;i++) FIC >> w[i];

    D = new double * [N];
    for(i=0;i<N;i++) D[i] = new double [N];

    DT = new double * [N];
    for(i=0;i<N;i++) DT[i] = new double [N];     

    int c = 0;
    for(x1 = 0; x1 < N ; x1 ++)
    for(x2 = 0; x2 < N ; x2 ++)
     {
		FIC >> d;				
        D[x1][x2] = d; 
        DT[x1][x2] = 2.0*d;
        c++;
     }

    FIC.close();
    /*
    for(i=0;i<N;i++)
     { 
      printf("\n");
      for(j=0;j<N;j++) printf("%d ", int (D[i][j])); 
     }
     */
 }
void AssignMemery()
{
    int i,j;  
	
    p = new int [N];
    bestp = new int [N];
    SizeG = new double [K];
	TabuTenure = new int *[N];
	for(i=0;i<N;i++) TabuTenure[i] = new int [K]; 
	
	PairSet = new Pair [N*K];
	
    Delta_Matrix = new double *[N];
       for(i=0;i<N;i++) Delta_Matrix[i] = new double [K]; 
    Delta_Matrix1 = new double *[N];
       for(i=0;i<N;i++) Delta_Matrix1[i] = new double [K];    
       
    Delta = new double *[N];
       for(i=0;i<N;i++) Delta[i] = new double [K];  
       
    Delta1 = new double *[N];
       for(i=0;i<N;i++) Delta1[i] = new double [K];      
     
	CS.p = new int [N];
	NS.p = new int [N];
	GS.p = new int [N]; 
	OS.p = new int [N]; 
    
	CS.SizeG = new double [K];
	NS.SizeG = new double [K];
	GS.SizeG = new double [K];
	OS.SizeG = new double [K];
	
	Neighbors = new Neighborhood [N*(N-1)/2 + N*K ]; 
}

void ReleaseMemery()
{    
     int i;
	 
     delete [] p; p = NULL;
	 delete [] bestp; bestp = NULL; 
	 delete [] SizeG; SizeG = NULL; 
	 
	 delete [] CS.p; CS.p = NULL;
	 delete [] CS.SizeG; CS.SizeG = NULL; 
	 delete [] GS.p; GS.p = NULL;
	 delete [] GS.SizeG; GS.SizeG = NULL; 
	 delete [] NS.p; NS.p = NULL;
	 delete [] NS.SizeG; NS.SizeG = NULL; 
	 delete [] OS.p; OS.p = NULL;
	 delete [] OS.SizeG; OS.SizeG = NULL; 
	 
	 delete [] LB; LB = NULL; 
	 delete [] UB; UB = NULL; 
	 delete [] Neighbors; Neighbors = NULL; 
	 
	 for(i=0;i<N;i++)
	 { 
	   delete [] Delta_Matrix[i]; Delta_Matrix[i] = NULL ;
       delete [] Delta_Matrix1[i]; Delta_Matrix1[i] = NULL ;
       delete [] Delta[i]; Delta[i] = NULL ;
       delete [] Delta1[i]; Delta1[i] = NULL ;
	   delete [] D[i]; D[i] = NULL; 
       delete [] DT[i]; DT[i] = NULL; 
	 }
   
}
/******************************************************************************************/
/*********************************    3. OutPuting Results   ******************************/
/******************************************************************************************/
int Proof(Solution &S)
{   
    int i,j;
    double ff;
    int flag ;
    ff = 0.0;
    for( i = 0 ; i < N ; i++ )  
	   for( j = i+1 ; j < N ; j++ )
	      {
	          if(S.p[i] == S.p[j])
               {   	
				  ff += D[i][j];
	           }
	      }
    S.cost = ff; 
    for(i=0;i<K;i++) S.SizeG[i] = 0.0; 
    for(i=0;i<N;i++) S.SizeG[S.p[i]] +=  w[i]; 
    flag = 1;
    for(i=0; i < K; i++) if(S.SizeG[i] < LB[i] - 1.0e-10 || S.SizeG[i]> UB[i] + 1.0e-10) { flag = 0; break;} 
    return flag; 
}

void Outputing(Solution &S, char *filename)
{
    int i;int r;
	FILE *fp; 
	char buff[80];
	r= rand()%1000;
	if(Proof(S)==0) return;
    sprintf(buff,"%s.sol",filename);
    fp=fopen(buff,"a+");
    fprintf(fp,"N = %d  G = %d  f = %lf\n", N , K, S.cost); 
    for(i=0;i<K;i++)
    fprintf(fp,"%lf   %lf   %lf\n", LB[i], UB[i], S.SizeG[i]); 
    printf("\n");
    for(i=0;i<N;i++)
    fprintf(fp,"%5.4d   %5.3d\n",i, S.p[i]); 
	fclose(fp);
}

void Out_results(double best , double ave,  double worst, double AvgTime, double deviation, char *filename, char instance[])
{
    int i;
	FILE *fp; 
	char buff[80];
    sprintf(buff,"%s",filename);
    fp = fopen(buff,"a+");
    fprintf(fp,"%s   %lf   %lf   %lf   %lf   %lf\n", instance, best, ave, worst, deviation, AvgTime); 
	fclose(fp); 
}

double Deviation(double arr[],int n)
{
    int i;
    double sum = 0,tmp = 0, x_avg;
    for(i = 0; i < n; ++i) sum += arr[i];
    x_avg = sum / n;
    for(i = 0; i < n; ++i)  tmp += (arr[i] - x_avg)*(arr[i] - x_avg);
    return sqrt(tmp/n); 
} // standard deviation  
/******************************************************************************************/
/****************   4. Constructive Heuristics for Initial Solution   *********************/
/******************************************************************************************/
void RandomInitiaSol(int p[], double SizeG[])
{     
     int i,j;
     int p1;
     int count, c1, Nc;
     int *Flag = new int [N];
     double *SizeGroup = new double [K];
     int *G = new int [K];
     int *VN = new int [N]; 
     for(i=0;i<K;i++) SizeGroup[i] = 0.0;
     for(i=0;i<N;i++) Flag[i] = 0;
     
     count = 0; 
     Nc = 0;
     for(i=0;i<K;i++) if(SizeGroup[i] < LB[i]) G[count++] = i;  
     for(j=0;j<N;j++) if(Flag[j]==0) VN[Nc++] = j; 
     while(count > 0)
     {  
        while(1)
        { 
          p1 = VN[rand()% Nc]; 
          c1 = G[rand()%count]; 
          if(SizeGroup[c1] + w[p1] <= UB[c1])
          {
            p[p1] = c1; 
            Flag[p1] = 1; 
            SizeGroup[c1] += w[p1]; 
            break;
          } 
        }    
        count = 0; 
        Nc = 0;
        for(i=0;i<K;i++) if(SizeGroup[i] < LB[i]) G[count++] = i; 
        for(j=0;j<N;j++) if(Flag[j]==0) VN[Nc++] = j; 
     }
       
     count = 0; Nc = 0;
     for(i=0;i<K;i++) if(SizeGroup[i] < UB[i]) G[count++] = i;  
     for(j=0;j<N;j++) if(Flag[j]==0) VN[Nc++] = j; 
     while(Nc > 0)
     {
        while(1)
        {
          p1 = VN[rand()% Nc]; 
          c1 = G[rand()%count]; 
          if(SizeGroup[c1] + w[p1] <= UB[c1])
          {
            p[p1] = c1; 
            Flag[p1] = 1; 
            SizeGroup[c1] += w[p1]; 
            break;
          }  
        }   
        count = 0; Nc = 0; 
        for(i=0;i<K;i++) if(SizeGroup[i] < UB[i]) G[count++] = i;  
        for(j=0;j<N;j++) if(Flag[j]==0) VN[Nc++] = j; 
     } 

     for(i=0;i<K;i++) SizeG[i] = SizeGroup[i];
     delete [] SizeGroup; SizeGroup= NULL;
     delete [] Flag; Flag = NULL; 
     delete [] G; G = NULL;
     delete [] VN; VN = NULL; 
    // printf("finish construction \n");
    // for(i=0;i<K;i++)printf("%lf   %lf   %lf\n",LB[i], UB[i], SizeG[i]);
}

void RandomInitiaSol1(int p[], double SizeG[])
{     
     int i,j;
     int p1;
     int count, c1, Nc;
     int *Flag = new int [N];
     double v_max; int j_max;
     double *SizeGroup = new double [K];
     int *G = new int [K];
     int *VN = new int [N]; 
     for(i=0;i<K;i++) SizeGroup[i] = 0.0;
     for(i=0;i<N;i++) Flag[i] = 0;
     
     count = 0; 
     v_max = -99999.0;
     for(i=0;i<K;i++) if(SizeGroup[i] < LB[i]) G[count++] = i;  
     for(j=0;j<N;j++) if(Flag[j]==0 && w[j] > v_max) { v_max = w[j];j_max = j; }
     while(count > 0)
     {  
        while(1)
        { 
          p1 = j_max; 
          c1 = G[rand()%count]; 
          if(SizeGroup[c1] + w[p1] <= UB[c1])
          {
            p[p1] = c1; 
            Flag[p1] = 1; 
            SizeGroup[c1] += w[p1]; 
            break;
          } 
        }    
     count = 0; 
     v_max = -99999.0;
     for(i=0;i<K;i++) if(SizeGroup[i] < LB[i]) G[count++] = i;  
     for(j=0;j<N;j++) if(Flag[j]==0 && w[j] > v_max) { v_max = w[j];j_max = j; }
     }
   //  for(i=0;i<K;i++) printf("%lf  %lf ",LB[i], SizeGroup[i]);   
     count = 0; 
     v_max = -99999.0;
     Nc = 0; 
     for(i=0;i<K;i++) if(SizeGroup[i] < UB[i]) G[count++] = i;  
     for(j=0;j<N;j++) if(Flag[j]==0 && w[j] > v_max) { v_max = w[j];j_max = j; }
     for(j=0;j<N;j++) if(Flag[j]==0) Nc++;
     while(Nc > 0)
     {
        while(1)
        {
          p1 = j_max; 
          c1 = G[rand()%count]; 
          if(SizeGroup[c1] + w[p1] <= UB[c1])
          {
            p[p1] = c1; 
            Flag[p1] = 1; 
            SizeGroup[c1] += w[p1]; 
            break;
          } 
         // for(j=0;j<K;j++) printf("%lf  ",SizeGroup[j]);  printf("%d",count);
         // for(j=0;j<N;j++) if(Flag[j]==0) printf(" %lf ", w[j]); printf("\n");
        }   
      count = 0; 
      v_max = -99999.0;
      Nc = 0; 
      for(i=0;i<K;i++) if(SizeGroup[i] < UB[i]) G[count++] = i;  
      for(j=0;j<N;j++) if(Flag[j]==0 && w[j] > v_max) { v_max = w[j]; j_max = j; }
      for(j=0;j<N;j++) if(Flag[j]==0) Nc++;   
     } 

     for(i=0;i<K;i++) SizeG[i] = SizeGroup[i];
     delete [] SizeGroup; SizeGroup= NULL;
     delete [] Flag; Flag = NULL; 
     delete [] G; G = NULL;
     delete [] VN; VN = NULL; 
    // printf("\n finish construction \n");
    // for(i=0;i<K;i++)printf("%lf   %lf   %lf\n",LB[i], UB[i], SizeG[i]);
}


void GRASPInitiaSol(int p[],double SizeG[])
{    
   int i,j,k;
   int v,g;
   int r,Nc;
   double p_max;
   int count = 0; 
   double **IK;
   double *RCL;
   int *Flag = new int [N];
   double *SizeGroup = new double [K];
   int *G = new int [K];
   int *VN = new int [N]; 
   IK = new double *[N];
   for(i=0;i<N;i++) IK[i] = new double [K];
   RCL = new double [N];
   
   for(i=0;i<K;i++) SizeGroup[i] = 0.0;
   for(i=0;i<N;i++) Flag[i] = 0;
   //a. assign randomly K vertices to K different subsets
   k = 0;
   while(k < K)
   {
      r = rand()%N;
      if(Flag[r]==0)
      {
        p[r] = k;
        SizeGroup[k] += w[r];
        Flag[r] = 1;
        k ++;  
      }
   }
   // b. Assin the vertices to subsets such that the lower bound cosntraints are satisifed. 
   for(k=0;k<K;k++)
   { 
     
      while(SizeGroup[k] < LB[k])
      {  
        p_max = -9999999.0;                
        for(i=0;i<N;i++) RCL[i] = 0;
        for(i=0;i<N;i++)
          if(Flag[i]==0)
          { 
            for(j=0;j<N;j++)
             if(Flag[j] ==1 && p[j]==k) RCL[i] += D[i][j];
            if(RCL[i]/w[i]> p_max ) p_max = RCL[i]/w[i]; 
          }
          
        count = 0;   
        for(i=0;i<N;i++)
        {
          if(Flag[i]==0 && RCL[i]/w[i] >= 0.8*p_max)
          {
            VN[count++] = i; 
          }
        }
        
        r = VN[rand()%count];
        p[r] = k;
        SizeGroup[k] += w[r];
        Flag[r] = 1;  
      }
     
   }
  // for(k=0;k<K;k++) printf("%lf %lf\n",SizeGroup[k],LB[k]);
   //c. assign the remaining vertices
   for(i=0;i<N;i++)
      for(k=0;k<K;k++) IK[i][k] = 0.0; 
   for(i=0;i<N;i++)
     {
        if(Flag[i]==0)
        for(j=0;j<N;j++)
        if(Flag[j]==1)
        {
            IK[i][p[j]] += D[i][j];            
        }
     }
     
   Nc =0; 
   for(j=0;j<N;j++) if(Flag[j]==0) Nc++; 
   while(Nc>0)
   {
     p_max = -99999.0;
     for(i=0;i<N;i++)
      for(j=0;j<K;j++)
      {
         if( (IK[i][j]/w[i] > p_max) && (Flag[i]==0) && (SizeGroup[j] + w[i] <= UB[j])) p_max = IK[i][j]/w[i];
      }
     count = 0; 
     for(i=0;i<N;i++)
      for(j=0;j<K;j++)
      {
         if( (IK[i][j]/w[i] >= 0.8*p_max) && (Flag[i]==0) && (SizeGroup[j] + w[i] <= UB[j]) ) { PairSet[count].i = i; PairSet[count].g = j; count++;}
      }
     
     r = rand() % count; 
     v = PairSet[r].i; 
     g = PairSet[r].g;
     
     p[v] = g;
     SizeGroup[g] += w[v];
     Flag[v] = 1;
     for(j=0;j<N;j++)
     if(Flag[j]==0) IK[j][g] += D[j][v];  
     
     Nc =0; for(j=0;j<N;j++) if(Flag[j]==0) Nc++;       
   }
   for(k=0;k<K;k++) SizeG[k] = SizeGroup[k];
   
   delete [] RCL;
   for(i=0;i<N;i++) delete [] IK[i];
   delete [] IK; 
   delete [] SizeGroup; SizeGroup = NULL;
   delete [] Flag; Flag = NULL; 
   delete [] G; G = NULL;
   delete [] VN; VN = NULL; 
   
  // double sum1,sum2;
  // sum1= sum2 =0;
  // for(i=0;i<N;i++) sum1+= w[i];
  // for(j=0;j<K;j++) sum2+= SizeG[j];
  // printf("%lf   %lf\n",sum1,sum2);
}

void GRASPInitiaSol1(int p[],double SizeG[])
{    
   int i,j,k;
   int v,g;
   int r,Nc;
   double p_max;
   int count = 0; 
   double **IK;
   double *RCL;
   int *Flag = new int [N];
   double *SizeGroup = new double [K];
   int *G = new int [K];
   int *VN = new int [N]; 
   IK = new double *[N];
   for(i=0;i<N;i++) IK[i] = new double [K];
   RCL = new double [N];
   
   for(i=0;i<K;i++) SizeGroup[i] = 0.0;
   for(i=0;i<N;i++) Flag[i] = 0;
   //a. assign randomly K vertices to K different subsets
   k = 0;
   while(k < K)
   {
      r = rand()%N;
      if(Flag[r]==0)
      {
        p[r] = k;
        SizeGroup[k] += w[r];
        Flag[r] = 1;
        k ++;  
      }
   }
   // b. Assin the vertices to subsets such that the lower bound cosntraints are satisifed. 
   for(k=0;k<K;k++)
   { 
     
      while(SizeGroup[k] < LB[k])
      {  
        p_max = -9999999.0;                
        for(i=0;i<N;i++) RCL[i] = 0;
        for(i=0;i<N;i++)
          if(Flag[i]==0)
          { 
            for(j=0;j<N;j++)
             if(Flag[j] ==1 && p[j]==k) RCL[i] += D[i][j];
            if(RCL[i]> p_max ) p_max = RCL[i]; 
          }
          
        count = 0;   
        for(i=0;i<N;i++)
        {
          if(Flag[i]==0 && RCL[i] >= 0.6*p_max)
          {
            VN[count++] = i; 
          }
        }
        
        r = VN[rand()%count];
        p[r] = k;
        SizeGroup[k] += w[r];
        Flag[r] = 1;  
      }
     
   }
  // for(k=0;k<K;k++) printf("%lf %lf\n",SizeGroup[k],LB[k]);
   //c. assign the remaining vertices
   for(i=0;i<N;i++)
      for(k=0;k<K;k++) IK[i][k] = 0.0; 
   for(i=0;i<N;i++)
     {
        if(Flag[i]==0)
        for(j=0;j<N;j++)
        if(Flag[j]==1)
        {
            IK[i][p[j]] += D[i][j];            
        }
     }
     
   Nc =0; 
   for(j=0;j<N;j++) if(Flag[j]==0) Nc++; 
   while(Nc>0)
   {
     p_max = -99999.0;
     for(i=0;i<N;i++)
      for(j=0;j<K;j++)
      {
         if( (IK[i][j] > p_max) && (Flag[i]==0) && (SizeGroup[j] + w[i] <= UB[j])) p_max = IK[i][j];
      }
     count = 0; 
     for(i=0;i<N;i++)
      for(j=0;j<K;j++)
      {
         if( (IK[i][j] >= 0.6*p_max) && (Flag[i]==0) && (SizeGroup[j] + w[i] <= UB[j]) ) { PairSet[count].i = i; PairSet[count].g = j; count++;}
      }
     
     r = rand() % count; 
     v = PairSet[r].i; 
     g = PairSet[r].g;
     
     p[v] = g;
     SizeGroup[g] += w[v];
     Flag[v] = 1;
     for(j=0;j<N;j++)
     if(Flag[j]==0) IK[j][g] += D[j][v];  
     
     Nc =0; for(j=0;j<N;j++) if(Flag[j]==0) Nc++;       
   }
   for(k=0;k<K;k++) SizeG[k] = SizeGroup[k];
   
   delete [] RCL;
   for(i=0;i<N;i++) delete [] IK[i];
   delete [] IK; 
   delete [] SizeGroup; SizeGroup = NULL;
   delete [] Flag; Flag = NULL; 
   delete [] G; G = NULL;
   delete [] VN; VN = NULL; 
   
}
void GreedyInitiaSol(int p[],double SizeGroup[])
{    
    
     int i,j,v,g;
     int Nc;
     int r, g_max;
     int count, start_index, cur_index;
     int tot_number;
     int sum = 0;
     int *Flag;
     double *SumG, MaxSumG; 
     
     SumG = new double [K]; 
     for(g=0;g<K;g++) SizeGroup[g] = 0;
     Flag = new int [N];
     for(i=0;i<N;i++) Flag[i] = 0;
     
     //a. assign randomly K elements to K distinct groups
     for(g=0;g<K;g++)
     {
         while(1)
         {
             r = rand()%N;
             if(Flag[r]==0)
             {
                p[r] = g;
                Flag[r] = 1; 
                SizeGroup[g] += w[r];
                break; 
             }
         }
     }
  
     //b. assign greedily the elements to satisfy the lower bound constraints of groups 
     Nc = 0; 
     for(g = 0; g < K; g++) if(SizeGroup[g] < LB[g]) Nc++;
      
     while(Nc > 0)
     {  
        
        cur_index = rand()%N; 
        do
        {
           cur_index = (cur_index + 1)%N;           
        }while(Flag[cur_index]); 
        
        for(g=0;g<K;g++) SumG[g] = 0.0;
        
        for(j=0;j<N;j++)
           if(Flag[j]==1) SumG[p[j]] += D[cur_index][j];
        //for(g=0;g<K;g++) SumG[g] /=  w[cur_index];
        for(g=0;g<K;g++) SumG[g] /= (SizeGroup[g] + w[cur_index]);
        MaxSumG = -999999.0; 
        for(g=0; g<K; g++)
         if(SizeGroup[g] < LB[g] && SumG[g] > MaxSumG)
         {
            MaxSumG = SumG[g];
            g_max = g; 
         }
        
        p[cur_index] = g_max;
        Flag[cur_index] = 1;
        SizeGroup[g_max] += w[cur_index]; 
        
        Nc = 0; 
        for(g = 0; g < K; g++) if(SizeGroup[g] < LB[g]) Nc++;
     }
    // for(g=0;g<K;g++)printf("%lf   %lf\n",LB[g],SizeGroup[g]);
     //c. assign the remaining the elements 
     Nc = 0; 
     for(i = 0; i < N; i++) if(Flag[i] ==0) Nc++;
     while(Nc > 0)
     {
        cur_index = rand()%N; 
        do
        {
          cur_index = (cur_index + 1)%N;           
        }while(Flag[cur_index]); 
        
        for(g=0; g<K; g++) SumG[g] = 0.0; 
        
        for(j=0;j<N;j++)
           if(Flag[j]==1) SumG[p[j]] += D[cur_index][j];
        //for(g=0;g<K;g++) SumG[g] /=  w[cur_index];
        for(g=0;g<K;g++) SumG[g] /= (SizeGroup[g] + w[cur_index]);
        MaxSumG = -999999.0; 
        for(g=0; g<K; g++)
         if(SizeGroup[g] + w[cur_index] <= UB[g] && SumG[g] > MaxSumG)
         {
            MaxSumG = SumG[g];
            g_max = g; 
         } 
        
        p[cur_index] = g_max;
        Flag[cur_index] = 1;
        SizeGroup[g_max] += w[cur_index]; 
        
       Nc = 0; for(i = 0; i < N; i++) if(Flag[i] ==0) Nc++;
    }  
    
    delete [] Flag; Flag = NULL; 
    delete [] SumG; SumG = NULL;    
}
/******************************************************************************************/
/**********    5.   Variable Neighborhood Desecent method      **********************/
/******************************************************************************************/
void BuildNeighbors()
{
     int i,j,g;
     int count;
     int SN = N*(N-1)/2 + N*K;
     count = 0; 
     for(i=0;i<N;i++)
       for(g=0;g<K;g++)
       {
          Neighbors[count].type = 1;
          Neighbors[count].v = i ;
          Neighbors[count].g = g; 
          count ++; 
       } 
     for(i=0;i<N;i++)
       for(j=i+1;j<N;j++)
       {
          Neighbors[count].type = 2;
          Neighbors[count].x = i; 
          Neighbors[count].y = j;
          count ++; 
       }         
}

//2.1 Clear delta matrix
void Clear_Delta_Matrix( )
{
	int x, g ;
	f = 0.0;
	for( x = 0 ; x < N ; x++ )
		for( g = 0 ; g < K ; g++ )
		Delta_Matrix[ x ][ g ] = 0.0 ; 
    for( x = 0 ; x < N ; x++ )
      for( g = 0 ; g < K ; g++ )
   	    TabuTenure[ x ][ g ] = 0 ; 
	return ;		
}

//2.2 Build delta matrix
void Build_Delta_Matrix()
{
	int i,j;
	Clear_Delta_Matrix( );
	for(i = 0; i < N ; i++ )
	   for( j = 0 ; j < N ; j++ )
	      Delta_Matrix[ i ][ p[j] ]  +=  D[i][j];		 
//	for(i = 0; i < N ; i++ )
//	   for( j = 0 ; j < K ; j++ )
//       Delta[i][j] = Delta_Matrix[ i ][ j ] - Delta_Matrix[ i ][ p[i] ];     
    f = 0.0;
    for( i = 0 ; i < N ; i++ )
      f += Delta_Matrix[i][p[i]];
    f = f/2;                  
  //printf("\n f= %lf ********** \n", f);
    return;   
}
//2.2 Update one move delta matrix
inline void One_Move_Update_Delta_Matrix1(int i, int g0, int g1)
{
    int x,j,k;
    
    for(j=0;j<N;j++)
     {  
           if(j!=i)
           {
             Delta_Matrix[ j ][ g0 ] -= D[i][j];
             Delta_Matrix[ j ][ g1 ] += D[i][j];
           } 
     }                                                         
	return ;     
}
void LocalSearch_N1(int p1[], double SizeG[], double *cost)
{
     int v, g, counter, old_g;
     double delt, ff; 
     
     counter = 0; 
     ff = *cost;
     while(counter < N*K)
     {
       for(v=0;v<N;v++)
       {
         for(g=0;g<K;g++)
           if( ( p1[v] != g ) && ( SizeG[ p1[v] ] - w[v] >= LB[ p1[v] ] ) && (SizeG[g] + w[v] <= UB[g] )  )
             {
               delt = Delta_Matrix[ v ][ g ] - Delta_Matrix[ v ][ p1[ v ] ]; 
               //printf("delt =%lf\n",delt);
               if(delt > 1.0e-11)
               {  
                  old_g = p1[v] ;              
				  One_Move_Update_Delta_Matrix1(v, old_g, g);
				  SizeG[old_g] = SizeG[old_g] - w[v];
                  SizeG[g] = SizeG[g] + w[v];
				  p1[v] = g;	 
			  	  ff += delt; 
			  	  counter = 0; 
			  	 // printf(" delt = %lf improve !\n",delt);
               }
               else counter++; 
             }
             else counter ++;
        if(counter >= N*K) break;   
       }
     }
     *cost = ff; 
     
    // printf("N1 f=%lf\n",ff); 
    // scanf("%d",&g);
}

int LocalSearch_N2_1(int p1[], double SizeG[], double *cost)
{
     int x, y, g, counter;
     double delt,ff; 
     int improve_cutoff; 
     int old_g, old_g1, swap;
     
     ff = *cost; 
     counter = 0; 
     improve_cutoff = 0;
     while(counter < N*(N-1)/2)
     {
        for(x=0;x<N;x++)
        {
           for(y=x+1;y<N;y++)
             if((p1[x] != p1[y])&&(SizeG[p1[x]] + w[y] - w[x] >= LB[p1[x]]) && (SizeG[p1[x]] + w[y] - w[x] <= UB[p1[x]]) && (SizeG[p1[y]] + w[x] - w[y] >= LB[p1[y]])&&(SizeG[p1[y]] + w[x] - w[y] <= UB[p1[y]]) )
              {  
              // delt = Delta[ x ][ p1[y] ] + Delta[ y ][ p1[x] ] - DT[x][y];   
               delt =  (Delta_Matrix[ x ][ p1[y] ] - Delta_Matrix[ x ][ p1[ x ] ]) + (Delta_Matrix[ y ][ p1[x] ] - Delta_Matrix[ y ][ p1[ y ] ]) - DT[x][y];
               if(delt > 1.0e-10)
               {
                old_g  = p1[ x ]; 
                old_g1 = p1[ y ];   
                SizeG[old_g]  += (w[y]-w[x]);
                SizeG[old_g1] += (w[x]-w[y]);     
				One_Move_Update_Delta_Matrix1( x, old_g, old_g1);
				
				swap    = p1[ x ];
				p1[ x ] = p1[ y ];
				One_Move_Update_Delta_Matrix1( y, old_g1, old_g);
				p1[ y ] = swap; 
                  
				ff += delt; 
			   //printf("swap delt = %lf\n",delt);
				counter = 0; 
				improve_cutoff ++; 
				if(improve_cutoff >= 1) {*cost = ff; return 1;}
               }
               else counter ++;
             }
             else counter ++;
          if(counter >= N*(N-1)/2) break; 
        }
     }
    *cost = ff; 
   // printf("N2 f=%lf\n",ff);
    if(improve_cutoff > 0) return 1; 
    else return 0; 
}

void VND2(int partition[], double SizeG[], double *cost)
{ 
     int i,rez;
     double fc;
     for(i=0;i<N;i++) p[i] = partition[i];
	
     for(i=0;i<K;i++) SizeG[i] = 0.0;
	 for(i=0;i<N;i++) SizeG[p[i]] += w[i];
	 
     Build_Delta_Matrix(); 
     fc = f;
     do
     {
        LocalSearch_N1(partition, SizeG, &fc);
        rez =  LocalSearch_N2_1(partition, SizeG, &fc);
     }while(rez);     
     *cost = fc;  
}

/******************************************************************************************/
/***********************************    6. Initial solution    ****************************/
/******************************************************************************************/
void InitialSol(Solution &S)
{
     
     GRASPInitiaSol1(S.p, S.SizeG); 
    
}

/******************************************************************************************/
/*****************************      7. GRASP Method     **********************************/
/******************************************************************************************/
void GRASP()
{
    int i,j;
    int L; 
    double f_c;
    starting_time = clock();
    GS.cost = -9999999; 
  
    while(1.0*(clock()- starting_time)/CLOCKS_PER_SEC  <  Time_limit)
    {  
      InitialSol(CS); 
      VND2(CS.p, CS.SizeG, &CS.cost); 
      if(CS.cost > GS.cost)
      {
        for(i=0;i<N;i++) GS.p[i] = CS.p[i]; 
        for(j=0;j<K;j++) GS.SizeG[j] = CS.SizeG[j]; 
        GS.cost = CS.cost;  
        total_time = ( clock() - starting_time )/CLOCKS_PER_SEC ; 
       // printf("GRASP f = %lf  time=%lf s \n", GS.cost,  total_time); 
      }
     
    }    	
}

/******************************************************************************************/
/********************    8. Main Function for Multiple Runs     **************************/
/******************************************************************************************/
int main(int argc, char *argv[])
{
    int i,j;
    int i1,j1;
    int seed; 
    double Sum;
    const int  Times = 20; 
    double F[Times], FM[Times];
    double Ctime[Times];
    double AvgTime = 0.0;
    double F_best= -99999999, F_worst = 999999999, F_ave = 0.0, deviation, deviation1;
    double FM_best= 99999999, FM_worst = -999999999, FM_ave = 0.0;
    seed = time(NULL) % 1000000 ;
	srand( seed );	
	
	 
   // File_Name = "RanReal240_01.txt";
   // Output_File_Name = "ss.txt";
   // Solution_File = "RanReal240_01.txt.sol";
     
     
    File_Name = argv[1]; 
    //Solution_File = argv[2];
    Output_File_Name = argv[2];
	//Time_limit = atoi(argv[4]); 
    Output_File_Name = "GRASP.txt";
    Solution_File = "FMS.txt"; 
    
    inputing();  
    AssignMemery();
    Time_limit = 1.0*N;
    BuildNeighbors();
    
    Sum = 0.0; 
    for(i=0; i<N; i++)
      for(j = i+1; j<N; j++)
      Sum += D[i][j]; 
    
    OS.cost = -9999999.0; 
    for(j=0;j<Times;j++) {F[j] = 0.0; Ctime[j] = 0; }
    for(i=0; i < Times; i++)
    { 
      GRASP();
      if(Proof(GS)) 
      { 
        F[i] = GS.cost; 
        FM[i] = 2.0*(Sum - F[i]);
        Ctime[i] = total_time; 
        if(F[i]> OS.cost)
        {
            for(i1=0;i1<N;i1++) OS.p[i1] = GS.p[i1]; 
            for(j1=0;j1<K;j1++) OS.SizeG[j1] = GS.SizeG[j1]; 
            OS.cost = GS.cost;  
        }
      }
     // printf("f[%d] = %lf \n",i, F[i]);
    }
    for(i=0;i<Times;i++)   
    {
       if(F[i] > F_best )  F_best = F[i];    
       if(F[i] < F_worst)  F_worst = F[i];  
       F_ave += F[i];   
       AvgTime += Ctime[i];
    }
    for(i=0;i<Times;i++)   
    {
       if(FM[i] < FM_best )  FM_best = FM[i];    
       if(FM[i] > FM_worst)  FM_worst = FM[i];  
       FM_ave += FM[i];   
    }
    F_ave   /=  Times; 
    FM_ave   /=  Times; 
    AvgTime /=  Times; 
    deviation = Deviation(F,Times);
    //deviation1 = Deviation(FM,Times);
    Out_results(F_best, F_ave, F_worst, AvgTime, deviation, Output_File_Name, File_Name); 
    //Out_results(FM_best, FM_ave, FM_worst, AvgTime, deviation1, Solution_File, File_Name);
    Outputing(OS, File_Name);
    ReleaseMemery();
    return 0;
}
