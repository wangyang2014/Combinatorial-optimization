#include <cstdlib>
#include <iostream>
using namespace std;
// FAP_PR.cpp : Defines the entry point for the console application.

//#pragma once
#define WIN32_LEAN_AND_MEAN		// 从 Windows 头中排除极少使用的资料
#include <stdio.h>
//#include <tchar.h>
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
//#include <conio.h>
#include <ctype.h>

/******************************************************/
/********************* 1.  Global variables ***********/
/******************************************************/
#define number_pop  30
#define number_elit  1
#define max_iteration 5000

typedef struct Adjacent_Matrix{
	int neighbor ;
	struct Adjacent_Matrix *next ;
}Adjacent ;

typedef struct POP_Class{
	int *p;
	int value;
	int distance;
	double G_x;
}POP_Class;
typedef struct PairSet{
	int i;
	int j;
	int number ; 
}Pair; 

typedef struct Element{
	int delta;
	int number ; 
}Element;

Adjacent * *A_Matrix ;         
Adjacent *p1, *q1; 
char * File_Name;
char Output_File_Name[30];
int N, K, G_K;  // node number and color number
unsigned int f, f_best ;
double total_time, starting_time; 
int * Color; // color array for each vertex
int * Best_Color ; // color array for each vertex
int * Best_Color_so_far;
int ** Delta_Matrix;  // incremental matrix for each vertex being colored each color
int ** Delta_Matrix1;
int ** Edge;   // adjacent matrix
int ** D;  //distance matrix
unsigned int ** P;   //constraint matrix
int ** TabuTenure;  // tabu tenure
int ** Freq; //
int maxFreq; 
int * V; //vertex 
int *Neibor_number;
POP_Class pop[number_pop];
POP_Class solution_best; 
POP_Class off_spring;
POP_Class pop_elit[number_elit];
POP_Class pop_ref[ number_pop - number_elit];
int **Pair_Set; 
Pair *pair_s; 
float *value_best, *value_ave;
const int alpha = 1000;
const double gema = 0.33; 
const int times =10; 
float *t1;
int obj_value[times];
int obj_best;
int obj_worst;
double obj_ave; 
double Tot_time;
double time_ave;

/***************************************************/
/******************** 2. Inputing date *************/
/***************************************************/

void inputing()
{
	int i, x, y, x1, x2,d,p;
	ifstream FIC;
	FIC.open(File_Name);
	if ( FIC.fail() )
	{
		cout << "### Erreur open, File_Name " << File_Name << endl;
		exit(0);
	}
	char StrReading[100];
	FIC >> StrReading;
	if ( FIC.eof() )
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}
	//int nb_vtx=0 ;
	int nb_edg=-1, max_edg=0;
	while ( ! FIC.eof() )
	{
		char bidon[50];
		if ( strcmp(StrReading, "p" )==0 )
		{
			FIC >> bidon >> N >> nb_edg;
			//cout << "Number of vertexes = " << N << endl;
			//cout << "Number of edges = " << nb_edg << endl;
            
			Color = new int[N];
			Best_Color = new int[N];      
			Best_Color_so_far= new int [N];
			for( x = 0 ; x < N ; x++ ) 
			{
				Color[x] = rand() % K ;
				// Move_Freq[ x ] = 0 ;
			} 
			
			Edge=new int*[N];
			for (x = 0 ; x < N ; x++ ) 
				Edge[x]= new int[N];
			D=new int*[N];
			for (x = 0 ; x < N ; x++ ) 
				D[x]= new int[N];
			P=new unsigned int*[N];
			for (x = 0 ; x < N ; x++ ) 
				P[x]= new unsigned int[N];
			
			A_Matrix = new Adjacent *[N];
			for( i = 0 ; i < N; i ++ )
			{
				A_Matrix[ i ] = new Adjacent ;
				A_Matrix[ i ]->neighbor = 0 ; 
				A_Matrix[ i ]->next = NULL ;
			}  
			
			Delta_Matrix=new int*[N];
			for (x=0;x<N;x++) Delta_Matrix[x]=new int[K];
			Delta_Matrix1=new int*[N];
			for (x=0;x<N;x++) Delta_Matrix1[x]=new int[K];

			TabuTenure=new int*[N];
			for (x=0;x<N;x++) TabuTenure[x]=new int[K];
			
			Freq=new int*[N];
			for (x = 0; x < N; x++) Freq[x]=new int[K];
			
			for (x=0;x<N;x++)
                for (y=0;y<N;y++)
				{
					Edge[x][y] = 0;
					D[x][y] = -1; 
					P[x][y] = 0;
				} 
		}
		
		if ( strcmp(StrReading, "e")==0)
		{
			FIC >> x1 >> x2 >>d >> p ;
			// cout << x1 <<"  "<< x2 <<"  "<<d<<" "<<p<<  endl;
			//x1--; x2--;  // one shift
			if ( x1<0 || x2<0 || x1>=N || x2 >=N )
			{
				cout << "### Error of node : x1="
					<< x1 << ", x2=" << x2 << endl;
				exit(0);
			}
			if(x1!=x2)
			{
				Edge[x1][x2]=Edge[x2][x1] = 1;
				D[x1][x2] = D[x2][x1] = d;
				P[x1][x2] = P[x2][x1] = ((p > 100000) ? (100000) : (p) );
				max_edg++;
				// add x2 to x1's neighbor list
				p1 = A_Matrix[ x1 ] ;
				A_Matrix[ x1 ]->neighbor ++ ;
				while( p1->next != NULL )
					p1 = p1->next;   
				q1 = new Adjacent ;
				q1->neighbor = x2 ;
				q1->next = NULL ;  
				p1->next = q1 ;
				
				// add x1 to x2's neighbor list
				p1 = A_Matrix[ x2 ] ;
				A_Matrix[ x2 ]->neighbor ++ ;
				while( p1->next != NULL )
                    p1 = p1->next;          
				q1 = new Adjacent ;
				q1->neighbor = x1 ;
				q1->next = NULL ;  
				p1->next = q1;
			}
		}
		
		
		FIC >> StrReading;
     }
	 
	 
	 //for(i=0;i<N;i++) printf(" %d\n ",A_Matrix[ i ]->neighbor);
    // cout << "Density = " << (float) 2 * max_edg/(N*(N-1)) << endl;
     if ( 0 && max_edg != nb_edg )
     {
		 cout << "### Error de lecture du graphe, nbre aretes : annonce="
			 << nb_edg << ", lu=" << max_edg  << endl;
		 exit(0);
     }
	 
     FIC.close();
}

void Assign_Memery(POP_Class pop[number_pop],POP_Class *solution_best, POP_Class *off_spring) 
{
	int i,j; 
    for(i=0;i<number_pop; i++)
		pop[i].p= new int [N]; 
	for(i=0;i<number_elit; i++)
		pop_elit[i].p= new int [N]; 
	for(i=0;i<number_pop-number_elit; i++)
		pop_ref[i].p= new int [N]; 
	(*solution_best).p = new int [N]; 
	(*off_spring).p = new int [N];  
    Pair_Set = new int * [number_pop]; 
	for(i=0;i<number_pop;i++) 
		Pair_Set[i] = new int [number_pop]; 
	
    for(i=0; i<number_pop; i++)
		for(j=0;j<number_pop;j++) Pair_Set[i][j] = 0; 
		pair_s =new Pair [number_pop*number_pop/2]; 
		value_best = new float [max_iteration];
		for(i=0;i<max_iteration;i++) value_best[i]=0.0;
		value_ave = new float [max_iteration];
		for(i=0;i<max_iteration;i++) value_ave[i]=0.0;
		t1 = new float [max_iteration];
		for(i=0;i<max_iteration;i++) t1[i]=0.0;
		
}

void DeleteMemery(POP_Class pop[number_pop],POP_Class *solution_best, POP_Class *off_spring) 
{
    int i; 
    for(i=0;i<number_pop; i++)
		delete [] pop[i].p; 
    for(i=0;i<number_elit; i++)
		delete [] pop_elit[i].p;
    for(i=0;i<number_pop-number_elit; i++)
		delete [] pop_ref[i].p;
	delete [] (*solution_best).p; 
	delete [] (*off_spring).p ;  
    for( i=0; i<number_pop; i++ ) delete [] Pair_Set[i]; 
	delete [] pair_s; 
	delete [] value_best;
	delete [] value_ave; 
	delete [] t1; 
}

/**************************************************/
/*****************3. Outputing*******************/
/**************************************************/

void WriteDate(int *C, char *filename, int K,int f)
{ 
	int i;
	FILE *fp; 
	char buff[80];
    sprintf(buff,"%s -------- %d.txt",filename, K);
    fp=fopen(buff,"a+");
    fprintf(fp,"\n K=%d     f=%d\n ",K,f);
	for(i=0;i<N;i++)
		fprintf(fp,"%d ",C[i]);	
	fclose(fp);
}//end writedate function

void WriteDate1(float *C,char *filename,int q)
{ 
	int i;
	FILE *fp; 
	char buff[80];
    sprintf(buff,"%s - %lf - %d.txt",filename, gema, alpha);
    fp=fopen(buff,"a+");
    fprintf(fp,"%d   %d \n", q, alpha);
	for(i=0;i<3000;i++)
		fprintf(fp,"%d      %lf \n", i, 1.0*C[i]/times);	
	fclose(fp); 
}//end writedate function
void WriteDate2(float *C,char *filename,int q)
{ 
	int i;
	FILE *fp; 
	char buff[80];
    sprintf(buff,"%s -- %d.txt",filename, alpha);
    fp=fopen(buff,"a+");
    fprintf(fp,"%d   %d \n", q, alpha);
	for(i=0;i<3000;i++)
		fprintf(fp,"%d     %lf \n", i, 1.0*C[i]/times);	
	fclose(fp); 
}//end writedate function
void Output_Results(int succ)
{
	
    FILE *fp ;
	// Total_Time = (double)(clock() - Starting_Time)/ 1000 ;
    fp = fopen(Output_File_Name, "a+"); 
    fprintf(fp," K=%d  number_succ= %d  total_time=%lf   average_time=%lf \n", K, succ,total_time,total_time/succ);
    fclose(fp) ;
    return ;
}
void Out_Put(int f_best, int f_worst, double f_ave,  double time_ave, double tot_time,char *filename )
{ 
	int i;
	FILE *fp; 
	char buff[80];
    sprintf(buff,"%s--%d.txt",filename, K);
    fp=fopen(buff,"a+");
    fprintf(fp,"K = %d\n ",K);
    fprintf(fp,"f_best = %d\n ",f_best);
    fprintf(fp,"f_worst = %d\n ",f_worst);
    fprintf(fp,"f_ave = %f\n ",f_ave);
    fprintf(fp,"time_ave  = %f\n ",time_ave);
    fprintf(fp,"tot_time  = %f\n ",tot_time);
	fclose(fp);
}//end writedate function
/*********************************************/
/*************** 4. Tabu Search **************/
/*********************************************/

void Clear_Delta_Matrix( )
{
	int x, v ;
	f = 0;
	for( x = 0 ; x < N ; x++ )
		for( v = 0 ; v < K ; v++ )
			Delta_Matrix[ x ][ v ] = 0 ;
		
		for( x = 0 ; x < N ; x++ )
			for( v = 0 ; v < K ; v++ )
			{
				TabuTenure[ x ][ v ] = 0 ;
				Freq[ x ][ v ] = 0;
			}
			maxFreq = 1;
			return ;     
}

//2.2 Build delta matrix
void Build_Delta_Matrix( )
{
	int i, j,s ;
	Clear_Delta_Matrix( ) ;
	
	for( i = 0 ; i < N ; i++ )
		for(s=0; s<K; s++)
			for( j = 0 ; j <N ; j++ )
			{
				if( (i!=j) && Edge[ i ][ j ] != 0 )
				{      
					if(abs(Color[ j ] - s)<=D[i][j])   Delta_Matrix[ i ][ s ]  +=  P[i][j];
				}
			}
			f=0;
			for( i = 0 ; i < N ; i++ )
				for( j = i+1 ; j <N ; j++ )
				{
					if((Edge[ i ][ j ] == 1)&&abs(Color[ i ] - Color[ j ])<=D[i][j])
					{   
						
						f += P[i][j];
					}
				}
		//		printf("f = %d ********** \n", f);
				return ;        
}

void One_Move_Update_Delta_Matrix(int i, int v0, int v1)
{
	int j , s;int start,end;
	p1 = A_Matrix[ i ] ;
	while( p1->next != NULL )    
	{ 
		p1 = p1->next ;
		j = p1->neighbor ;
		if((v0-D[i][j])<0)start=0;
		else start=v0-D[i][j];
		if((v0+D[i][j])>K-1)end=K-1;
		else end= v0+D[i][j];
		for(s=start; s<=end; s++)
		{  
			Delta_Matrix[ j ][ s ] -= P[i][j];
		}
		if((v1-D[i][j])<0)start=0;
		else start=v1-D[i][j];
		if((v1+D[i][j])>K-1)end=K-1;
		else end= v1+D[i][j];
		for(s=start; s<=end; s++)
		{
			Delta_Matrix[ j ][ s ] += P[i][j];
		}
	}                                                                        
	return ;     
}
void One_Move_Update_Delta_Matrix1(int i, int v0, int v1)
{
	int j , s;int start,end;
	p1 = A_Matrix[ i ] ;
	while( p1->next != NULL )    
	{ 
		p1 = p1->next ;
		j = p1->neighbor ;
		if((v0-D[i][j])<0)start=0;
		else start=v0-D[i][j];
		if((v0+D[i][j])>K-1)end=K-1;
		else end= v0+D[i][j];
		for(s=start; s<=end; s++)
		{  
			Delta_Matrix1[ j ][ s ] -= P[i][j];
		}
		if((v1-D[i][j])<0)start=0;
		else start=v1-D[i][j];
		if((v1+D[i][j])>K-1)end=K-1;
		else end= v1+D[i][j];
		for(s=start; s<=end; s++)
		{
			Delta_Matrix1[ j ][ s ] += P[i][j];
		}
	}                                                                        
	return ;     
}
int One_Move_Tabu_Search(int Color_into[], int *value)
{
	int i ;
	int non_improve = 0 ;  // the stop condition of TS
	int num_tabu_best, num_best ;  // the number of tabu neighbors and non-tabu neighbors
	int best_x[ 50 ], best_v[ 50 ], x,  v;
	int tabu_best_x[ 50 ], tabu_best_v[ 50 ];
	int iter;
	int tabu_best_delta, best_delta, delt ;
	int old_color ;
	int select ;
	int num ;
	int R;
	
	double total_time, starting_time;
	starting_time = clock();
	for(i=0;i<N;i++)
		Color[i] = Color_into[i];
	
	
	Build_Delta_Matrix( );
	f_best = f ;
//	printf("\n\n");
//	cout << endl << "One_Move :       iter       f       f_best      time " << endl;
//	cout << "---------------------------------------------------" << endl; 
//	R = 20 + rand()%60 +1; 
//	R = int (0.3*N);
//	R = 50; 
	iter = 0 ; 
	while( non_improve < 1000 )
	{
		num = 0 ;
		tabu_best_delta = 9999999 ; 
		best_delta = 9999999 ;
		num_tabu_best = 0 ; 
		num_best = 0 ;	
		for( x = 0 ; x < N ; x++ )
            if( Delta_Matrix[ x ][ Color[ x ] ] )
			{
                for( v = 0 ; v < K ; v++ )
					if( v != Color[ x ] )
                    {
						num ++ ;
						delt = Delta_Matrix[ x ][ v ] - Delta_Matrix[ x ][ Color[ x ] ];
						if( TabuTenure[ x ][ v ] <= iter ) // if this is not tabued 
                        {
							if( delt < best_delta )
							{
								best_x[ 0 ] = x ; 
								best_v[ 0 ] = v ;
								best_delta = delt ; 
								num_best = 1 ;
							}
							else if( delt == best_delta && num_best < 50 )
							{
								best_x[ num_best ] = x ; 
								best_v[ num_best ] = v ;
								num_best ++ ;
							}
						}                                                    
						else if( TabuTenure[ x ][ v ] > iter )// if it is tabu 
						{ 
							if( delt < tabu_best_delta  )
							{
								tabu_best_x[ 0 ] = x ; 
								tabu_best_v[ 0 ] = v ; 
								tabu_best_delta = delt ; 
								num_tabu_best = 1 ;
							}
							else if( delt == tabu_best_delta && num_tabu_best < 50 )
							{
								tabu_best_x[ num_tabu_best ] = x ; 
								tabu_best_v[ num_tabu_best ] = v ; 
								num_tabu_best ++ ;
							}                               
						}
					}
			}        
			//choose the tabu best move if the tab aspiration criterion is satisfied
			if( ( num_tabu_best > 0 && tabu_best_delta < best_delta && ( f + tabu_best_delta < f_best ) ) || num_best == 0 )  // aspiration criterion 
			{
				f += tabu_best_delta ;
				select = rand( ) % num_tabu_best ;  
				// to select the best tabu move according to the second level criteria
				
				old_color = Color[ tabu_best_x[ select ] ] ;              
				One_Move_Update_Delta_Matrix( tabu_best_x[ select ], old_color, tabu_best_v[ select ] );
				Color[ tabu_best_x[ select ] ] = tabu_best_v[ select ] ;
				
				Freq[tabu_best_x[ select ]][old_color]+=1;
				if( Freq[tabu_best_x[ select ]][old_color]>maxFreq) maxFreq=Freq[tabu_best_x[ select ]][old_color];
				
				TabuTenure[ tabu_best_x[ select ] ][ old_color ] =  10+ rand()%3+ 1 +(int)(40*(1.0*Freq[tabu_best_x[ select ]][old_color ]/maxFreq)); // + (int)co_freq * Move_Freq[ tabu_best_x[ select ] ]; 
				TabuTenure[ tabu_best_x[ select ] ][ old_color ] += iter ; 
			} 
			else //choose the best non tabu move
			{
				f += best_delta ; 
				select = rand( ) % num_best ;    
				old_color = Color[ best_x[ select ] ] ;              
				One_Move_Update_Delta_Matrix( best_x[ select ], old_color, best_v[ select ] );
				Color[ best_x[ select ] ] = best_v[ select ] ;
				
				Freq[best_x[ select ]][old_color]+=1;
			    if( Freq[best_x[ select ]][old_color]>maxFreq) maxFreq=Freq[best_x[ select ]][old_color];
				
				TabuTenure[ best_x[ select ] ][ old_color ] =  10 + rand()%3 + 1 +(int)(40*(1.0*Freq[best_x[ select ]][old_color ]/maxFreq)); // + (int) co_freq * Move_Freq[ best_x[ select ] ] ; // (int)( alpha * f ) ;
				TabuTenure[ best_x[ select ] ][ old_color ] += iter ; 
			} 
			
			iter ++ ;
			//Total_Iterations ++ ;
			
			total_time = (clock() - starting_time )/CLOCKS_PER_SEC ;
			
			if( f <= f_best )
			{
				if( f < f_best )
				{
					f_best = f;
					for( i = 0 ; i < N ; i ++ )
						Best_Color[ i ] = Color[ i ] ;
					
	//				printf("\n One_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
					
					non_improve = 0 ;
				}  
				else if ( f == f_best )   
					non_improve ++ ;
				if( f_best == 0 )
				{
					//printf("\n == %d      %d      %d  ", iter, f, f_best);
					for(i=0;i<N;i++)
						Color_into[i] = Best_Color[ i ];
					*value = f_best; 
					return f_best;
				}  
			}  
			else   
				non_improve ++ ;   
     }
	 for(i=0;i<N;i++)
		 Color_into[i] = Best_Color[ i ];
	 *value = f_best; 
     return f_best;
}
int One_Move_Desent(int Color_into[], int *value)
{
	int i ;
	int non_improve = 0 ;  // the stop condition of TS
	int num_best ;  // the number of tabu neighbors and non-tabu neighbors
	int best_x[ 50 ], best_v[ 50 ], x,  v;
    // int tabu_best_x[ 50 ], tabu_best_v[ 50 ];
	int iter;
	int  best_delta, delt ;
	int old_color ;
	int select ;
	int num ;
	
	// double total_time, starting_time;
	// starting_time = clock();
	for(i=0;i<N;i++)
		Color[i] = Color_into[i];
	
	
	Build_Delta_Matrix( );
	f_best = f ;
    // printf("\n\n");
    // cout << endl << "One_Move :       iter       f       f_best      time " << endl;
	//  cout << "---------------------------------------------------" << endl; 
	
	iter = 0 ; 
	best_delta = -1;
	while( best_delta < 0 )
	{
		num = 0 ;
		
		best_delta = 9999999 ;
		
		num_best = 0 ;	
		for( x = 0 ; x < N ; x++ )
            if( Delta_Matrix[ x ][ Color[ x ] ] )
			{
                for( v = 0 ; v < K ; v++ )
					if( v != Color[ x ] )
                    {
						num ++ ;
						delt = Delta_Matrix[ x ][ v ] - Delta_Matrix[ x ][ Color[ x ] ];
						
						if( delt < best_delta )
						{
							best_x[ 0 ] = x ; 
							best_v[ 0 ] = v ;
							best_delta = delt ; 
							num_best = 1 ;
						}
						else if( delt == best_delta && num_best < 50 )
						{
							best_x[ num_best ] = x ; 
							best_v[ num_best ] = v ;
							num_best ++ ;
						}
						
						
					}
			}        
			
			if(best_delta < 0)
			{
				f += best_delta ; 
				select = rand( ) % num_best ;    
				old_color = Color[ best_x[ select ] ] ;              
				One_Move_Update_Delta_Matrix( best_x[ select ], old_color, best_v[ select ] );
				Color[ best_x[ select ] ] = best_v[ select ] ;
			} 
			else break;
			
			iter ++ ;
			//Total_Iterations ++ ;
			
			// total_time = (clock() - starting_time )/CLOCKS_PER_SEC ;
			
			if( f <= f_best )
			{
				if( f < f_best )
				{
					f_best = f;
					for( i = 0 ; i < N ; i ++ )
						Best_Color[ i ] = Color[ i ] ;
					
					// printf("\n One_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
					
					
				}  
				
				if( f_best == 0 )
				{
					//printf("\n == %d      %d      %d  ", iter, f, f_best);
					for(i=0;i<N;i++)
						Color_into[i] = Best_Color[ i ];
					*value = f_best; 
					return f_best;
				}  
			}  
			
     }
	 for(i=0;i<N;i++)
		 Color_into[i] = Best_Color[ i ];
	 *value = f_best; 
     return f_best;
}

int One_Move_Steepest_Desent(int Color_into[], int *value)
{
	int i ;
	int  num_best ;  // the number of best neighbors.
	int best_x[ 50 ], best_v[ 50 ], x, v;
	int iter;
	int best_delta, delt ;
	int old_color;
	int select ;
    
	
	//	 double total_time, starting_time;
	// starting_time = clock();
	for(i=0;i<N;i++)
		Color[i] = Color_into[i];
	Build_Delta_Matrix( );
	f_best = f ;
    // printf("\n\n");
//	cout << endl << "One_Move :       iter       f       f_best      time " << endl;
//	cout << "---------------------------------------------------" << endl; 
	
	iter = 0 ; 
	best_delta = -1;
	while(best_delta<0)
	{        
		best_delta = 9999999;
		num_best = 0 ;	
		for( x = 0 ; x < N ; x++ )
            if( Delta_Matrix[ x ][ Color[ x ] ] )
			{
                for( v = 0 ; v < K ; v++ )
					if( v != Color[ x ] )
					{
						delt = Delta_Matrix[ x ][ v ] - Delta_Matrix[ x ][ Color[ x ] ];
						if(delt < best_delta )
                        {
							best_x[ 0 ] = x ; 
							best_v[ 0 ] = v ;
							best_delta = delt ; 
							num_best = 1 ;
						}
						else if( delt == best_delta && num_best < 50 )
						{
							best_x[ num_best ] = x ; 
							best_v[ num_best ] = v ;
							num_best ++ ;
						}                                                                 
					}
			}  
			if(best_delta<0)
			{
				f += best_delta ;
				select = rand( ) % num_best ;  
				// to select the best  move according to the  steepest decesent criteria.   
				old_color = Color[ best_x[ select ] ] ;              
				One_Move_Update_Delta_Matrix(best_x[ select ], old_color, best_v[ select ] );
				Color[ best_x[ select ] ] = best_v[ select ] ;
			}
			
			iter ++ ;
			//  total_time = (clock() - starting_time )/CLOCKS_PER_SEC ;
			if( f <= f_best )
			{
				if( f < f_best )
				{
					f_best = f;
					for( i = 0 ; i < N ; i ++ )
						Best_Color[ i ] = Color[ i ] ;
					//printf("\n One_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
				}  
				if( f_best == 0 )
				{
					//  printf("\n == %d      %d      %d  ", iter, f, f_best);
					for(i=0;i<N;i++)
						Color_into[i] = Best_Color[ i ];
					*value = f_best; 
					return f_best;
				}  
			}  
			
     }
	 for(i=0;i<N;i++)
		 Color_into[i] = Best_Color[ i ];
	 *value = f_best; 
     return f_best;
}


/*********************************************************/
/**************** 5. Path_Relinking operator ************/
/********************************************************/

void Path_Relinking(int *x,int *y, int *off_spring)
{
    double gema = 0.4;
	int f_min = 9999999;
	int i;
	int j; 
	int *NC;
	int *PV;
	int *FI;
	int k=0;
	int r=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(i=0;i<N;i++) NC[i]=0; 
	PV = new int [N+1];
	for(i=0;i<=N;i++)PV[i]=-1;
	FI = new int [N+1];
	for(i=0;i<=N;i++)FI[i]=0;
	for(i=0;i<N;i++)Color[i]=x[i];
	Build_Delta_Matrix( );
	FI[0]=f;
    for(i=0;i<N;i++) if(x[i]!=y[i]){ r++; NC[i]=1;}
	for(i=1;i <= r;i++)
	{
		min_delta = 99999999;
		for(j=0;j<N;j++) 
		{    
			if(NC[j]!=0)
			{
				delta = Delta_Matrix[j][y[j]] - Delta_Matrix[j][x[j]]; 
				if(delta < min_delta) { k = j; min_delta = delta;}
			}	  
		}
		FI[i] = FI[i-1] + min_delta;
		One_Move_Update_Delta_Matrix( k, x[k], y[k]); 
		NC[k] = 0; 
		PV[i] = k;
		
	}
	int s = 0; 
	f_min=9999999;
	
	for(i=(int)(gema*r);i<(int)((1.0-gema)*r);i++) 
	{   
		if(FI[i] < f_min){ f_min = FI[i]; s = i;}
	}
//	for(i=1;i<=r;i++)   printf ("f(x[%d])= %d\n ", i, FI[i]);
  //  printf("f_min = %d \ n ", f_min); 
	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=s;i++) off_spring[PV[i]]=y[PV[i]]; 
	
	delete [] NC;
	delete [] FI;
	delete [] PV;
}
void Mixed_Path_Relinking(int *x,int *y, int *off_spring)
{
	
	int f_min = 9999999;
	int count1 = 0, count2 = 0; 
	int i;
	int j; 
	int *NC;
	int *PV;
	int *PV1;
	int *FI;
	int *FI1;
	int k=0;
	int r=0;int r1=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(i=0;i<N;i++) NC[i]=0; 
	PV = new int [N+1];
	for(i=0;i<=N;i++)PV[i]=-1;
	PV1 = new int [N+1];
	for(i=0;i<=N;i++)PV1[i]=-1;
	FI = new int [N+1];
	for(i=0;i<=N;i++)FI[i]=0;
	FI1 = new int [N+1];
	for(i=0;i<=N;i++)FI1[i]=0;

	for(i=0;i<N;i++)Color[i]=x[i];
	Build_Delta_Matrix( );
	for(i=0;i<N;i++)
		for(j=0;j<K;j++) 
			Delta_Matrix1[i][j] = Delta_Matrix[i][j]; 
	FI1[0]=f;
	

	for(i=0;i<N;i++)Color[i]=y[i];
	Build_Delta_Matrix( );
	FI[0]=f;
	
	for(i=0;i<N;i++) if(x[i]!=y[i]){ r++; NC[i]=1;}
	for(i=1;i <= r;i++)
	{
		min_delta = 99999999;

		if(i%2==1)
		{
			for(j=0;j<N;j++) 
			{    
				if(NC[j]!=0)
				{
					delta = Delta_Matrix1[j][y[j]] - Delta_Matrix1[j][x[j]]; 
					if(delta < min_delta) { k = j; min_delta = delta;}
				}	  
			}
			FI1[count1+1] = FI1[count1] + min_delta;
			One_Move_Update_Delta_Matrix1( k, x[k], y[k]); 
			NC[k] = 0; 
			PV1[count1] = k;
			count1++;
			
		}
		else 
		{
			for(j=0;j<N;j++) 
			{    
				if(NC[j]!=0)
				{
					delta = Delta_Matrix[j][x[j]] - Delta_Matrix[j][y[j]]; 
					if(delta < min_delta) { k = j; min_delta = delta;}
				}	  
			}
			FI[count2 + 1] = FI[count2] + min_delta;
			One_Move_Update_Delta_Matrix( k, y[k], x[k]); 
			NC[k] = 0; 
			PV[count2] = k;
			count2++; 
			
		}

	}
	
//	for(i=0;i<count1;i++) printf("x[%d]= %d \n",i,FI1[i]);
	
//	for(i=count2;i>=0;i--) printf("y[%d]= %d \n",i,FI[i]);


	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=count1;i++) off_spring[PV1[i]] = y[PV1[i]]; 

//  printf("r=%d,   count1= %d \n", r, count1);
//	r1=0;
//	for(i=0;i<N;i++)if(off_spring[i] != x[i]) r1++;
//	printf("distance to x is %d \n",r1);
//	r1=0;
//	for(i=0;i<N;i++)if(off_spring[i] != y[i])r1++;
//	printf("distance to y is %d \n",r1);


	delete [] NC; NC=NULL;
	delete [] FI; FI=NULL;
	delete [] FI1; FI1= NULL;
	delete [] PV; PV = NULL;
	delete [] PV1; PV1= NULL; 
}
void Path_Relinking1(int *x,int *y, int *off_spring)
{
	// double gema = 0.35;
	int f_min = 9999999;
	int i;
	//	int j; 
	int *NC;
	int *PV;
	int *FI;
	int k=0;
	int r=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(i=0;i<N;i++) NC[i]=0; 
	PV = new int [N+1];
	for(i=0;i<=N;i++)PV[i]=-1;
	FI = new int [N+1];
	for(i=0;i<=N;i++)FI[i]=0;
	for(i=0;i<N;i++)Color[i]=x[i];
	Build_Delta_Matrix( );
	FI[0]=f;
    for(i=0;i<N;i++) if(x[i]!=y[i]){ r++; NC[i]=1;}
	for(i=1;i <= r;i++)
	{
		
		while(1)
		{   
			k = rand() % N;
			if(NC[k]!=0)
			{ 
				delta = Delta_Matrix[k][y[k]] - Delta_Matrix[k][x[k]];
				FI[i] = FI[i-1] + delta;
				One_Move_Update_Delta_Matrix( k, x[k], y[k]); 
				NC[k] = 0; 
				PV[i] = k; 
				break;
			}
			
		}
		
	}
	
	int s = 0; 
	f_min=9999999;
	
	for(i=(int)(gema*r);i<(int)((1.0-gema)*r);i++) 
	{   
		if(FI[i] < f_min){ f_min = FI[i]; s = i;}
	}
//	for(i=1;i<=r;i++) { if(i%10==0)printf("\n");  printf ("f(x[%d])= %d ", i, FI[i]);}
 //   printf("f_min = %d \ n ", f_min); 
	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=s;i++) off_spring[PV[i]]=y[PV[i]]; 
	
	delete [] NC;
	delete [] FI;
	delete [] PV;
}
void Path_Relinking2(int *x,int *y, int *off_spring)
{
	// double gema = 0.35;
	int f_min = 9999999;
	int i,j;
	//	int j; 
	int *NC;
	int *NC1;
	int *PV;
	int *FI;
	int k=0;
	int r=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(j=0;j<N;j++) NC[j]=0; 
	NC1 = new int [N+1];
    for(j=0;j<=N;j++) NC1[j]=0; 
	int count_number = 0; 
	PV = new int [N+1];
	for(j=0;j<=N;j++)PV[j]=-1;
	FI = new int [N+1];
	for(j=0;j<=N;j++)FI[j]=0;
	for(j=0;j<N;j++)Color[j]=x[j];
	Build_Delta_Matrix( );
	FI[0]=f;
    for(j=0;j<N;j++) if(x[j]!=y[j]){ r++; NC[j]=1;}
	for(i=1;i <= r;i++)
	{    
		count_number = 0;
		for(j=0;j<N;j++) if( NC[j]!=0 ) { NC1[count_number] = j; count_number++; }
		k = NC1[rand() % count_number]; 
		if(NC[k]!=0)
		{ 
			delta = Delta_Matrix[k][y[k]] - Delta_Matrix[k][x[k]];
			FI[i] = FI[i-1] + delta;
			One_Move_Update_Delta_Matrix( k, x[k], y[k]); 
			NC[k] = 0; 
			PV[i] = k; 	
		}
	//	printf("f= %d \n",FI[i]);
	}
	
	int s = 0; 
	f_min=9999999;
	
	for(i=(int)(gema*r);i<(int)((1.0-gema)*r);i++) 
	{   
		if(FI[i] < f_min){ f_min = FI[i]; s = i;}
	}

	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=s;i++) off_spring[PV[i]]=y[PV[i]]; 
	
	delete [] NC;
	delete [] FI;
	delete [] PV;
}
void Path_Relinking3(int *x,int *y, int *off_spring)
{
	// double gema = 0.35;
	int f_min = 9999999;
	int i,j;
	//	int j; 
	int *NC;
	int *NC1;
	int *PV;
	int *FI;
	int k=0;
	int r=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(j=0;j<N;j++) NC[j]=0; 
	NC1 = new int [N+1];
    for(j=0;j<=N;j++) NC1[j]=0; 
	int count_number = 0; 
	PV = new int [N+1];
	for(j=0;j<=N;j++)PV[j]=-1;
	FI = new int [N+1];
	for(j=0;j<=N;j++)FI[j]=0;
	for(j=0;j<N;j++)Color[j]=x[j];
	Build_Delta_Matrix( );
	FI[0]=f;
    for(j=0;j<N;j++) if(x[j]!=y[j]){ r++; NC[j]=1;}
	for(i=1;i <= r;i++)
	{    
		count_number = 0;
		for(j=0;j<N;j++) if( NC[j]!=0 ) { NC1[count_number] = j; count_number++; }
		k = NC1[rand() % count_number]; 
		if(NC[k]!=0)
		{ 
			delta = Delta_Matrix[k][y[k]] - Delta_Matrix[k][x[k]];
			FI[i] = FI[i-1] + delta;
			One_Move_Update_Delta_Matrix( k, x[k], y[k]); 
			NC[k] = 0; 
			PV[i] = k; 	
		}
		
	}
	
	int s = 0; 
	f_min=9999999;
	
	for(i=(int)(gema*r);i<(int)((1.0-gema)*r);i++) 
	{   
		if(FI[i] < f_min){ f_min = FI[i]; s = i;}
	}

	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=s;i++) off_spring[PV[i]]=y[PV[i]]; 
	
	delete [] NC;
	delete [] FI;
	delete [] PV;
}
void Randomized_Path_Relinking(int *x,int *y, int *off_spring)
{
	
	int f_min = 9999999;
	int i,j;
	int q;
	Element * CSL;
	CSL = new Element [N];
	for(i=0;i<N;i++) {CSL[i].delta =0;CSL[i].number =0;}
	int *NC;
	int *NC1;
	int *PV;
	int *FI;
	int k1=0;
	int r=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(j=0;j<N;j++) NC[j]=0; 
	NC1 = new int [N+1];
	for(j=0;j<=N;j++) NC1[j]=0; 
	int count_number = 0; 
	PV = new int [N+1];
	for(j=0;j<=N;j++)PV[j]=-1;
	FI = new int [N+1];
	for(j=0;j<=N;j++)FI[j]=0;
	for(j=0;j<N;j++)Color[j]=x[j];
	Build_Delta_Matrix( );
	FI[0]=f;
	for(j=0;j<N;j++) if(x[j]!=y[j]){ r++; NC[j]=1;}
	for(q=1;q <= r;q++)
	{    

	   // build a CSL
		count_number = 0;
		for(i=0;i<N;i++)
		{
			if(NC[i]!=0)
			{ 
				delta = Delta_Matrix[i][y[i]] - Delta_Matrix[i][x[i]];
				if(count_number< 10 )
				{
					CSL[count_number].delta = delta;
					CSL[count_number].number = i;
					count_number++;
				}
				else 
				{
					int kk=0;
					int delta_max=0;
					for(j=0;j<count_number;j++)
					{
						if(CSL[j].delta > delta_max){ delta_max = CSL[j].delta; kk = j; }
					}
					if(delta < delta_max) { CSL[kk].delta= delta; CSL[kk].number=i;}
					if(delta==delta_max&&(rand()%2==0)) { CSL[kk].delta= delta; CSL[kk].number=i;}
				}
			}

		}

		// Select a move from the CSL

		 k1 = rand()%count_number; 
		{ 
			FI[q] = FI[q-1] + CSL[k1].delta;
			One_Move_Update_Delta_Matrix(  CSL[k1].number, x[ CSL[k1].number], y[ CSL[k1].number]); 
			NC[CSL[k1].number] = 0; 
			PV[q] = CSL[k1].number; 	
		}

	//	printf("f = %d \n ",FI[q]);
	}

	int s = 0; 
	f_min=9999999;

	for(i=(int)(gema*r);i<(int)((1.0-gema)*r);i++) 
	{   
		if(FI[i] < f_min){ f_min = FI[i]; s = i;}
	}

	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=s;i++) off_spring[PV[i]]=y[PV[i]]; 
	//printf("r=%d\n",r);
	//printf("s=%d\n",s);
	//s=0;
	//for(i=0;i<N;i++) if(off_spring[i]!=x[i])s++;
	//printf("distance to X = %d\n",s);

	delete [] NC;
	delete [] FI;
	delete [] PV;
	delete [] CSL;
}

void Combination_operator_one(int *pop_x1, int *pop_x2, int *off_spring)
{
	int i;
	for(i=0; i<N; i++)
	{ 
		if(rand()%2==0) off_spring[i] = pop_x1[i];
		else off_spring[i] = pop_x2[i];
	}
}
void Mixed_Randomized_Path_Relinking(int *x,int *y, int *off_spring)
{

	int f_min = 9999999;
	int i,j;
	int q;
	int count1=0, count2=0; 
	Element * CSL;
	CSL = new Element [N];
	for(i=0;i<N;i++) { CSL[i].delta =0; CSL[i].number =0; }
	int *NC;
	int *NC1;
	int *PV;
	int *PV1;
	int *FI;
	int *FI1;
	int k=0;
	int r=0;
	int delta=0;
	int min_delta=99999999;
	NC = new int [N+1];
	for(j=0;j<N;j++) NC[j]=0; 
	
	int count_number = 0; 
	PV = new int [N+1];
	for(j=0;j<=N;j++)PV[j]=-1;
	PV1= new int [N+1];
	for(i=0;i<=N; i++) PV1[i]=0;
	FI = new int [N+1];
	for(j=0;j<=N;j++)FI[j]=0;
	FI1= new int [N+1];
	for(j=0;j<=N;j++)FI1[j]=0;

	for(i=0;i<N;i++)Color[i]=x[i];
	Build_Delta_Matrix( );
	for(i=0;i<N;i++)
		for(j=0;j<K;j++) 
			Delta_Matrix1[i][j] = Delta_Matrix[i][j]; 
	FI1[0]=f;

	for(i=0;i<N;i++)Color[i]=y[i];
	Build_Delta_Matrix( );
	FI[0]=f;

	for(j=0;j<N;j++) if(x[j]!=y[j]){ r++; NC[j]=1;}

	for(q=1;q <= r;q++)
	{    
		min_delta = 99999999;
		if(q % 2 == 1)
		{
			// build a CSL
			count_number = 0;
			for(i=0;i<N;i++)
			{
				if(NC[i]!=0)
				{ 
					delta = Delta_Matrix1[i][y[i]] - Delta_Matrix1[i][x[i]];
					if(count_number< 5)
					{
						CSL[count_number].delta = delta;
						CSL[count_number].number = i;
						count_number++;
					}
					else 
					{
						int kk=0;
						int delta_max=0;
						for(j=0;j<count_number;j++)
						{
							if(CSL[j].delta > delta_max){ delta_max = CSL[j].delta; kk = j; }
						}
						if(delta < delta_max) { CSL[kk].delta= delta; CSL[kk].number=i;}
						if(delta==delta_max&&(rand()%2==0)) { CSL[kk].delta= delta; CSL[kk].number=i;}
					}
				}

			}

			// Select a move from the CSL
			k =rand()%count_number; 
			{ 
				FI1[ count1 + 1 ] = FI1[ count1 ] + CSL[k].delta;
				One_Move_Update_Delta_Matrix1(  CSL[k].number, x[ CSL[k].number], y[ CSL[k].number] ); 
				NC[CSL[k].number] = 0; 
				PV1[count1] = CSL[k].number; 	
				count1++;
			}
		}
		else
		{
			count_number = 0;
			for(i=0;i<N;i++)
			{
				if(NC[i]!=0)
				{ 
					delta = Delta_Matrix[i][x[i]] - Delta_Matrix[i][y[i]];
					if(count_number < 5 )
					{
						CSL[count_number].delta = delta;
						CSL[count_number].number = i;
						count_number++;
					}
					else 
					{
						int kk=0;
						int delta_max=0;
						for(j=0;j<count_number;j++)
						{
							if(CSL[j].delta > delta_max){ delta_max = CSL[j].delta; kk = j; }
						}
						if(delta < delta_max) { CSL[kk].delta= delta; CSL[kk].number=i;}
						if(delta==delta_max&&(rand()%2==0)) { CSL[kk].delta= delta; CSL[kk].number=i;}
					}
				}

			}

			// Select a move from the CSL
			k =rand()%count_number; 
			{ 
				FI[ count2 + 1 ] = FI[count2] + CSL[k].delta;
				One_Move_Update_Delta_Matrix(  CSL[k].number, y[ CSL[k].number], x[ CSL[k].number]); 
				NC[CSL[k].number] = 0; 
				PV[count2] = CSL[k].number; 	
				count2++;
			}
		}
		
	}
	//for(i=0;i<count1;i++) printf(" f1(%d)  =  %d \n ",i, FI1[i]);
	//for(i=0;i<count2;i++) printf(" f2(%d)  =  %d \n ",i, FI[i]);

	for(i=0;i<N;i++) off_spring[i]=x[i];
	for(i=1;i<=count1;i++) off_spring[PV1[i]]=y[PV1[i]]; 


	delete [] NC;
	delete [] FI1;
	delete [] FI;
	delete [] PV;
	delete [] PV1; 
	delete [] CSL; 
}



/************************************************************/
/************* 6.population initilization*******************/
/***********************************************************/

void  pop_initilization() 
{
	
	int i,j;  
	int count;
	int f_max;
	int f_min=99999999;
	int k;
	int k_min;
	int s;
	POP_Class xx; 
	xx.p = new int [N];
	
	for(i=0;i<number_pop;i++)
		for(j=0;j<N;j++)   pop[i].p[j] = rand()%K;	
		for(i=0;i<number_pop;i++) 
		{ 
			One_Move_Tabu_Search(pop[i].p, &(pop[i].value)); 
			if(pop[i].value==0) {for(j=0;j<N;j++)solution_best.p[j] = pop[i].p[j]; solution_best.value = pop[i].value; delete [] xx.p; xx.p = NULL; return; }
			
		} 
		count=0;
		while(count < 2*number_pop)
		{   
			f_max = -1;
			k = 0 ;
			for(i=0;i<number_pop;i++) 
			{
				if(pop[i].value > f_max) { k=i; f_max = pop[i].value; }
			}
			for(j=0;j<N;j++)  xx.p[j] = rand()%K;
			One_Move_Tabu_Search(xx.p, &(xx.value)); 
			if(xx.value==0)
			{
				for(j=0;j<N;j++)solution_best.p[j] = xx.p[j];
				solution_best.value = xx.value; 
				delete [] xx.p;
				return; 
			}
			if(xx.value < f_max) 
			{
				for(j=0;j<N;j++)
					pop[k].p[j]=xx.p[j];
				pop[k].value = xx.value; 
			}
			
			count++; 
		}
		
		for(i=0;i<number_pop;i++) 
		{
			if(pop[i].value < f_min) { k_min = i; f_min = pop[i].value;  }
		}
		
		for(j=0;j<N;j++) solution_best.p[j] = pop[k_min].p[j];   
		solution_best.value = f_min; 
		
		
		for(i=0;i<number_pop;i++)
			for(j=i+1;j<number_pop;j++) Pair_Set[i][j] = 1;
			s=0; 
			for(i=0;i<number_pop;i++)
				for(j=i+1;j<number_pop;j++) 
					if(Pair_Set[i][j] ==1) {pair_s[s].i =i; pair_s[s].j =j; s++;} 
					
					pair_s[0].number = s;     
					
					
					for(j=0;j<N;j++)  off_spring.p[j] = rand()%K;  
					off_spring.value = 9999999; 
					
					delete [] xx.p; 
					
}

void  pop_initilization1() 
{
	
	int i,j;  
	int count;
	int f_max;
	int fp_max=-99999;
	int k;
	int s;
	int k_max;
	POP_Class xx; 
	xx.p = new int [N];
	
	for(i=0;i<number_pop;i++)
		for(j=0;j<N;j++)   pop[i].p[j] = rand()%K;	
		for(i=0;i<number_pop;i++) 
		{ 
			One_Move_Tabu_Search(pop[i].p, &(pop[i].value)); 
			if(pop[i].value==0) {for(j=0;j<N;j++) solution_best.p[j] = pop[i].p[j];solution_best.value = pop[i].value; delete [] xx.p;xx.p = NULL; return; }
			
		} 
		count=0;
		while(count < 2*number_pop)
		{   
			f_max = -1;
			k = 0 ;
			for(i=0;i<number_pop;i++) 
			{
				if(pop[i].value > f_max) { k=i; f_max = pop[i].value; }
			}
			for(j=0;j<N;j++)  xx.p[j] = rand()%K;
			One_Move_Tabu_Search(xx.p, &(xx.value)); 
			if(xx.value==0)
			{
				for(j=0;j<N;j++)solution_best.p[j] = xx.p[j];
				solution_best.value = xx.value; 
				delete [] xx.p; xx.p = NULL;
				return; 
			}
			if(xx.value < f_max) 
			{
				for(j=0;j<N;j++)
					pop[k].p[j]=xx.p[j];
				pop[k].value = xx.value; 
			}
			
			count++; 
		}
		
		for(i=0;i<number_pop;i++) 
		{
			if(pop[i].value > fp_max) { k_max = i; fp_max = pop[i].value;  }
		}
		
		for(j=0;j<N;j++) pop[k_max].p[j] = solution_best.p[j];   
		pop[k_max].value = solution_best.value ;
		
		
		for(i=0;i<number_pop;i++)
			for(j=i+1;j<number_pop;j++) Pair_Set[i][j] = 1;
			s=0; 
			for(i=0;i<number_pop;i++)
				for(j=i+1;j<number_pop;j++) 
					if(Pair_Set[i][j] ==1) {pair_s[s].i =i; pair_s[s].j =j; s++;} 
					
					pair_s[0].number = s;    
					
					for(j=0;j<N;j++)  off_spring.p[j] = rand()%K;  
					off_spring.value = 9999999; 
					
					delete [] xx.p; xx.p = NULL;
					
}


/************************************************************/
/********************* 7. Pool updating  *******************/
/***********************************************************/

int pop_updating(POP_Class pop[number_pop],POP_Class *off_spring,int *position)
{
	
	int i,j;
	int f_max = -999; 
	int k;
	int flag =0; 
	int count;

     
	
	for(i=0;i<number_pop;i++)
	{
		count =0;
		for(j=0;j<N;j++) if(pop[i].p[j]!=(*off_spring).p[j]) count++; 
		if( count < 0.2*N ) { return 0;} 	   
	}  
	
	for(i=0;i<number_pop;i++) 
	{
		if(pop[i].value > f_max) { k=i; f_max = pop[i].value; }
	}
	
	if((*off_spring).value < f_max) 
	{
		for(j=0;j<N;j++) pop[k].p[j]= (*off_spring).p[j]; 
		pop[k].value =  (*off_spring).value; 
		(*position) = k;
	}
	if((*off_spring).value < f_max) return 1;
	else return 0;  
}
int pop_updating1(POP_Class pop[number_pop],POP_Class *off_spring,int *position)
{
	
	int i,j;
	int f_max = -999999999; 
	int f_min = 999999;
	int k_max,k_min;
	int k_closest;
	int count;
	int similarity=99999990;
	
    for(i=0;i<number_pop;i++)
     {
       if(pop[i].value > f_max) { k_max = i; f_max = pop[i].value; } 
       if(pop[i].value < f_min) { k_min = i; f_min = pop[i].value; }            
     } 
	if((*off_spring).value > f_max) return 0; // don't need to update 
	
	for(i=0;i<number_pop;i++) 
	{
		count =0;
		for(j=0;j<N;j++) if(pop[i].p[j]!=(*off_spring).p[j]) count++; 
		if( count < similarity ) { similarity= count; k_closest=i;} 	   
	} // find out the most similar solution 
	
	if((*off_spring).value < f_min) 
	{
	    for(j=0;j<N;j++) pop[k_closest].p[j]= (*off_spring).p[j]; 
		pop[k_closest].value =  (*off_spring).value; 
		(*position) = k_closest;    
        return 1;              
    }// deplace the most simlilar solution by the offspring if the offspring is better than the best solution in the population
    if(((*off_spring).value < f_max)&&(similarity > 0.25*N))
    {
        for(j=0;j<N;j++) pop[k_closest].p[j]= (*off_spring).p[j]; 
		pop[k_closest].value =  (*off_spring).value; 
		(*position) = k_closest;    
        return 1;                                                
    }
    
   return 0;  
}
int pop_updating2(POP_Class pop[number_pop],POP_Class *off_spring,int *position)
{
	
	int i,j;
	int f_max = -99999; 
	int f_min = 99999;
	int k_max,k_min;
	int k_closest;
	int count;
	int similarity=99999;
	
    for(i=0;i<number_pop;i++)
     {
       if(pop[i].value > f_max) { k_max = i; f_max = pop[i].value; } 
       if(pop[i].value < f_min) { k_min = i; f_min = pop[i].value; }            
     } 
	if((*off_spring).value > f_max) return 0; // don't need to update 
	
	for(i=0;i<number_pop;i++) 
	{
		count =0;
		for(j=0;j<N;j++) if(pop[i].p[j]!=(*off_spring).p[j]) count++; 
		if( count < similarity ) { similarity= count; k_closest=i;} 	   
	} // find out the most similar solution 
	
	if((*off_spring).value < f_min) 
	{
	    for(j=0;j<N;j++) pop[k_closest].p[j]= (*off_spring).p[j]; 
		pop[k_closest].value =  (*off_spring).value; 
		(*position) = k_closest;    
        return 1;              
    }// deplace the most simlilar solution by the offspring if the offspring is better than the best solution in the population
    if(((*off_spring).value < pop[k_closest].value)&&(similarity <= 0.35*N))
    {
        for(j=0;j<N;j++) pop[k_closest].p[j]= (*off_spring).p[j]; 
		pop[k_closest].value =  (*off_spring).value; 
		(*position) = k_closest;    
        return 1;                                                
    }
    if(((*off_spring).value < f_max)&&(similarity > 0.35*N))
    {
        for(j=0;j<N;j++) pop[k_max].p[j]= (*off_spring).p[j]; 
		pop[k_max].value =  (*off_spring).value; 
		(*position) = k_max;    
        return 1;                                                
    }
    
   return 0;  
}
/***********************************************************/
/********************* 8. Updating PairSet *****************/
/**********************************************************/
void updating_PairSet(int x,int y)
{
	int i,j;
	int s=0;
	Pair_Set[x][y] = 0 ; 
	for(i=0;i<number_pop;i++)
		for(j=i+1;j<number_pop;j++) 
			if(Pair_Set[i][j] ==1) {pair_s[s].i =i; pair_s[s].j =j; s++;} 
			
			pair_s[0].number = s;  
}

void updating_PairSet_pop(int position)
{
	int i,j,s=0;
	for(i=position+1;i<number_pop;i++) Pair_Set[position][i] = 1;
	for(j=0;j<position;j++) Pair_Set[j][position] = 1; 
	
	for(i=0;i<number_pop;i++)
		for(j=i+1;j<number_pop;j++) 
			if(Pair_Set[i][j] ==1) {pair_s[s].i =i; pair_s[s].j =j; s++;} 
			
			pair_s[0].number = s;  
}
/***********************************************************/
/********************* 8.9 Average **********************/
/**********************************************************/
float AveragePOP()
{
	int i ;
	float f_ave=0;
	for(i=0;i<number_pop;i++) f_ave += 1.0*pop[i].value; 
    return (f_ave/number_pop); 
}

/***********************************************************/
/********************* 9. Rebuild population ************/
/**********************************************************/
void pop_elit_prev(POP_Class pop[number_pop],POP_Class pop_elit[number_elit])
{
 int i,j,k;
 POP_Class xx;
 xx.p = new int [N];
 for(i=0;i<number_pop;i++)
   for(j=i+1;j<number_pop;j++) 
   { 
      if(pop[i].value>pop[j].value) 
      {
          for(k=0;k<N;k++) xx.p[k] =  pop[i].p[k] ; 
          for(k=0;k<N;k++) pop[i].p[k] =  pop[j].p[k] ;    
          for(k=0;k<N;k++) pop[j].p[k] =  xx.p[k] ;      
          xx.value= pop[i].value;
          pop[i].value =  pop[j].value;
          pop[j].value =  xx.value;                                        
      }                                                    
   }
   
 for(i=0;i<number_elit;i++)
  {
    for(k=0;k<N;k++) pop_elit[i].p[k]= pop[i].p[k]; 
    pop_elit[i].value = pop[i].value; 
  }
   
 delete [] xx.p;    
 
}

void pop_ref_rebuild(POP_Class pop_ref[number_pop-number_elit])
{
    int i,j;  
	int count;
	int f_max;
	int k;
	POP_Class xx; 
	xx.p = new int [N];
	for( i=0; i < number_pop - number_elit ; i++ )
		for(j=0;j<N;j++)  pop_ref[i].p[j] = rand()%K;	
	for(i=0;i<number_pop-number_elit;i++) 
	{ 
	   One_Move_Tabu_Search(pop_ref[i].p, &(pop_ref[i].value)); 		
	} 
		count = 0;
		while(count < 3*number_pop )
		{   
			f_max = -1;
			k = 0 ;
			for(i=0;i < number_pop - number_elit;i++) 
			{
				if(pop_ref[i].value > f_max) { k=i; f_max = pop_ref[i].value; }
			}
			for(j=0;j<N;j++)  xx.p[j] = rand()%K;
			One_Move_Tabu_Search(xx.p, &(xx.value)); 
			if(xx.value < f_max) 
			{
				for(j=0;j<N;j++) pop_ref[k].p[j]=xx.p[j];
				pop_ref[k].value = xx.value; 
			}	
			count++; 
		}
	delete [] xx.p; xx.p = NULL;	
}

/***********************************************************/
/******************** 10. Main_Scheme  ******************/
/**********************************************************/

int  main(int argc, char **argv)
{
	int  seed ; 
	int x1,x2; 
    int f_min= 999999;
    int f_max= -99999;
    int k_min;
	int Iter;
	int IterationNumber = 0; 
	int s;
	int len;
	double START, END;
	int p; 
	int i,j;
	int position; 
	int number_rebuild =0;
	seed = time(NULL) % 1000000 ;
	srand( seed );	
	File_Name = argv[1]; 
	K = atoi(argv[2]); 
	//times = atoi(argv[3]); 	
	//File_Name = "45-17.ctr.txt"; 
	//printf("input k:\n");
	//scanf("%d", &K);	
	len = strlen(File_Name);
	strcpy(Output_File_Name, File_Name) ;
	Output_File_Name[len]='.';
	Output_File_Name[len+1]='t';
	Output_File_Name[len+2]='x';
	Output_File_Name[len+3]='t';
	Output_File_Name[len+4]='\0'; 
	inputing();	
    Assign_Memery(pop,&solution_best,&off_spring);     
	starting_time= clock();
	for(Iter=0; Iter< times; Iter++)
	{   
		IterationNumber = 0;
		number_rebuild = 0; 
		START = clock();
		for(i=0;i<number_pop;i++)
		for(j=0;j<N;j++)pop[i].p[j] = rand()%K;	
		for(i=0;i<number_pop;i++) 
		{ 
			One_Move_Tabu_Search(pop[i].p, &(pop[i].value));      		
		} //initilization for the first time
		IterationNumber = number_pop; 
		
		do{
		   
            pop_elit_prev(pop,pop_elit);	
            pop_ref_rebuild(pop_ref);
            
            IterationNumber = IterationNumber + 4*number_pop;
            
           for(i=0;i<number_elit;i++)
            {
                for(j=0;j<N;j++)pop[i].p[j]=pop_elit[i].p[j];
                pop[i].value =   pop_elit[i].value;                   
            }
            for(i=number_elit;i<number_pop;i++)
            {
                for(j=0;j<N;j++)pop[i].p[j] = pop_ref[i-number_elit].p[j];
                pop[i].value =   pop_ref[i-number_elit].value;                          
            }         
            //a. build the pairset table
           	for(i=0;i<number_pop;i++)
			for(j=i+1;j<number_pop;j++) Pair_Set[i][j] = 1;
			s=0; 
			for(i=0;i<number_pop;i++)
				for(j=i+1;j<number_pop;j++) 
					if(Pair_Set[i][j] ==1) {pair_s[s].i =i; pair_s[s].j =j; s++;} 
			pair_s[0].number = s;  
             //b. find the best solution 
             f_min=9999999;
             for(i=0;i<number_pop;i++)
              {
                if(pop[i].value<f_min){f_min= pop[i].value; k_min= i;}
              }
              for(j=0;j<N;j++)solution_best.p[j]=pop[k_min].p[j];
              solution_best.value= f_min; 
              //printf("best= %d \n",solution_best.value);     
              number_rebuild++;	
                   	
			//2. starting iterate  
			do{ 
			//	printf("size = %d\n",pair_s[0].number);
				p = rand()%pair_s[0].number; 
				x1= pair_s[p].i;
				x2= pair_s[p].j;
				updating_PairSet(x1,x2); 
				
				//1. x1------->x2
				Path_Relinking2(pop[x1].p,pop[x2].p,off_spring.p);				
				One_Move_Tabu_Search(off_spring.p, &(off_spring.value));  IterationNumber++;	  
				if(off_spring.value < solution_best.value) 
				{
                    for(j=0; j<N; j++)  solution_best.p[j] = off_spring.p[j];  
					solution_best.value = off_spring.value;  	
					//printf("inproved!******"); printf("best=%d\n",solution_best.value);
				}		
				if(pop_updating2(pop,&off_spring,&position)) 
				{  
					updating_PairSet_pop(position);
				}
				
				//2. x2------->x1
			    Path_Relinking2(pop[x2].p,pop[x1].p,off_spring.p);				
				One_Move_Tabu_Search(off_spring.p, &(off_spring.value)); IterationNumber++;	
				if(off_spring.value < solution_best.value) 
				{
                    for(j=0; j<N; j++)  solution_best.p[j] = off_spring.p[j];  
					solution_best.value = off_spring.value;  	
					//printf("*****inproved!******"); printf("best=%d \n",solution_best.value);
				}		
				if(pop_updating2(pop,&off_spring,&position)) 
				{  
					updating_PairSet_pop(position); 
				}	
                	      		  
			} while(pair_s[0].number!=0&&((clock()-START)/CLOCKS_PER_SEC) < 2400);   
	
          		  
		    //END = clock();
			  
		 } while (number_rebuild < 10000&&((clock()-START)/CLOCKS_PER_SEC) < 2400);  
		 obj_value[Iter] = solution_best.value;
		 WriteDate(solution_best.p,Output_File_Name, K, solution_best.value);
		 	 
 }//the number of runs.  
 
 Tot_time = (double)((clock()- starting_time)/CLOCKS_PER_SEC); 
 time_ave = Tot_time/times; 
 obj_best = 99999;
 obj_worst = -9999;
 obj_ave = 0.0; 
 for(i=0;i<times;i++)
 {    
       if(obj_value[i] > obj_worst)  obj_worst = obj_value[i];   
       if(obj_value[i] < obj_best )  obj_best  = obj_value[i];      
       obj_ave +=  (double)obj_value[i]; 
 }
 obj_ave /= times; 
 Out_Put(obj_best,obj_worst,obj_ave,time_ave,Tot_time,File_Name); 
 
 DeleteMemery(pop,&solution_best,&off_spring); 

// getchar(); 
 return 0; 
}




