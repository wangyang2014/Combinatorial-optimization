  
/*****************************************************************************\
 *                            Algorithm:  BITS-ECP                           *
 * This program demonstrates the use of the heuristic algorithm for solving  *
 * the equitable graph colouring problem (ECP). For more information about   *
 * this procedure (named as BITS-ECP), please email to: laixiangjing@gmail.com      *
\*****************************************************************************/ 
                                                                                                                                         
/*****************************************************************************/
/*****************    0. Header files and global varialbes       *************/
/*****************************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include<math.h>
#include<ctype.h>

using namespace std;

char * File_Name;
FILE *fp ;
char * outfilename;
int Seed ; 
int Total_Iterations ;
double starting_time, total_time, time_limit ;

typedef struct Adjacent_Matrix{
         int neighbor ;
         struct Adjacent_Matrix *next ;
         }Adjacent ;
typedef struct Solution{
        int *p;
        int *SizeG;
        int value;
        int k1;
        }Solution; 
typedef struct Statistics{
         int K_best;
         double Time_hit;
        }Statistics;
        
Adjacent * *A_Matrix ;         
Adjacent *p1, *q1; 

int N, K, G_K;  // node number and color number
int f, f_best ;
int Lbound;
int Ubound; 
int * Color; // color array for each vertex
int * SizeG;  
int * Best_Color ; // color array for each vertex
int ** Delta_Matrix;// incremental matrix for each vertex being colored each color
int ** Initial_Matrix; // for initialization of a k-coloring
int ** Edge;   // adjacent matrix
int ** TabuTenure;  // tabu tenure
int * Initial_Flag ;
Solution S, GS_best, SC;
Statistics Results[100]; // Statistics for computational results over multiple runs 
double Time_one_run_hit, Time_one_run_stating;
int T_max;
const int ALPHA = 100000; 

/*****************************************************************************/
/*****************          1. Initializing           ************************/
/*****************************************************************************/ 

//1.1 InPutting
void Initializing()
{
     int i, j, x, y, x1, x2; 
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
              Delta_Matrix = new int *[N];
              Initial_Matrix = new int *[N];         
              TabuTenure = new int *[N];
                 
              Edge=new int*[N];
              for (x = 0 ; x < N ; x++ ) 
                  Edge[x]= new int[N];
              
             
              A_Matrix = new Adjacent *[N];
              for( i = 0 ; i < N; i ++ )
                 {
                   A_Matrix[ i ] = new Adjacent ;
                   A_Matrix[ i ]->neighbor = 0 ; 
                   A_Matrix[ i ]->next = NULL ;
                 }  
              
              Initial_Flag = new int [N];
              for (x=0;x<N;x++)
                for (y=0;y<N;y++)
                    {
                      Edge[x][y]=0;
                    } 
   
			}
           
		   if ( strcmp(StrReading, "e")==0 )
           {
                 FIC >> x1 >> x2;
                 // cout << x1 << x2 << endl;
                 x1--; x2--;
                 if ( x1<0 || x2<0 || x1>=N || x2 >=N )
                 {
                       cout << "### Error of node : x1="
                            << x1 << ", x2=" << x2 << endl;
                       exit(0);
                 }
                 Edge[x1][x2]=Edge[x2][x1]=1;
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
                 p1->next = q1 ;
           }
           FIC >> StrReading;
     }
     //cout << "Density = " << (float) 2 * max_edg/(N*(N-1)) << endl;
     if ( 0 && max_edg != nb_edg )
     {
           cout << "### Error de lecture du graphe, nbre aretes : annonce="
                 << nb_edg << ", lu=" << max_edg  << endl;
           exit(0);
     }
   
     FIC.close();
}
void NewBound()
{     
      Lbound =  int (floor(1.0*N/K));
      Ubound =  int (ceil(1.0*N/K)) ;       
}

void free_memery()
{
     int i,j;
     for(i=0;i<N;i++) 
     {
      delete [] Delta_Matrix[i]; Delta_Matrix[i] = NULL;
      delete [] Initial_Matrix[i];Initial_Matrix[i]= NULL;
      delete [] TabuTenure[i]; TabuTenure[i] = NULL;
     } 
     delete [] SizeG; SizeG = NULL;
     delete [] S.p; S.p = NULL; 
     delete [] S.SizeG; S.SizeG = NULL;
     delete [] GS_best.p; GS_best.p = NULL;
     delete [] GS_best.SizeG; GS_best.SizeG = NULL;
}

void Initia_sol1(Solution &SS)
{
     int i,j,c1,c2,r;
     int *flag ;
     int *T ; 
     flag = new int [N]; 
     T = new int [N];
     for(i=0;i<N;i++) flag[i] = 0;  
     c2 = 0; 
     for(i=0;i<N;i++)
     {  
        c1=0;
        for(j=0;j<N;j++) if(flag[j]==0) { T[c1]=j; c1++;  }
        r = T[rand() % c1] ; 
        SS.p[r] = c2 ; 
        c2 = (c2 + 1) % K;  
        flag[r] = 1; 
     }
     for(i=0;i<K;i++) SS.SizeG[i]=0;
     for(i=0;i<N;i++) SS.SizeG[SS.p[i]]++; 
     //printf("**********************! \n");
     //for(i=0;i<K;i++) printf("%d ",SS.SizeG[i]); 
     //int sum =0; 
     //for(i=0;i<K;i++) sum += SS.SizeG[i];
     //printf("sum = %d \n",sum);
     
     delete [] flag; 
     delete [] T; 
}// To generate an initial solution for k-ECP

void Clear_Initial_Matrix()
{
	int x, v ;
	for(x=0; x < N; x++) Initial_Flag[x] = 0;
	for( x = 0 ; x < N ; x++ )
		for( v = 0 ; v < K ; v++ )
			Initial_Matrix[ x ][ v ] = 0 ;			
	return ;     
}

void Initia_sol2(int *patition, int *SizeG)
{
    int i,j,s,t,r,k;
    int flag; 
    int count; 
    int f;
    Adjacent *p;
    int *Inital_Set = new int [K];
    int CSL[50];
    Clear_Initial_Matrix();
    for(i=0;i<K;i++) Inital_Set[i]=-1;
    //a. choose randomly K vertices and assign them to different K classes.  
    for(s=0;s<K;s++)
    { 
     
      
      do{  
           flag = 1;
           r = rand()%N;  
		   t=0;
           while(t<s)
           {
               if(r==Inital_Set[t]) {flag=0;  break;}
			   t++;
           }
           if(flag == 1)  break;   
      }while(1);
      
     //printf("S[%d] = %d \n",s,r);
     
     Inital_Set[s] = r;
    }
    
	//for(k=0;k<K;k++)printf("%d    ",Inital_Set[k]);

    for(k=0;k<K;k++)
    {
       patition[Inital_Set[k]] = k; 
       Initial_Flag[Inital_Set[k]]=1;
       p = A_Matrix[Inital_Set[k]] ;
	  // printf("%d \n",Inital_Set[k]); 
       while( p->next != NULL )
       { 
           p = p->next; 
		   if(Initial_Flag[p->neighbor]==0) { Initial_Matrix[p->neighbor][k]++;}  //  printf("M[%d][%d] = %d \n",p->neighbor,k,Initial_Matrix[p->neighbor][k]); }
       }   
       //printf("****************************\n");            
    }
     
    //b. add the remaining vertices to classes in a one-by-one way
    int min;
    int nb;
    count = K; 
    while( count < N )
    {
        for( k = 0; k < K; k ++ )
        {
          min = 99999999;
          nb = 0;
		  if(count>=N) break;
          for(i=0;i<N;i++)
          {
            if(  Initial_Flag[i]==0  ) 
            {  
                if(Initial_Matrix[i][k] < min)
                { 
                    min = Initial_Matrix[i][k]; 
                    CSL[0] = i;
                    nb = 1; 
                }  
                if(Initial_Matrix[i][k] == min && nb < 30)
                {
                    CSL[nb]= i;
                    nb++;                
                }
              
            }
          }
  //     printf("min = %d \n",min);
         r = rand()% nb;
         patition[CSL[r]] = k; 
         Initial_Flag[CSL[r]] = 1; 
         p = A_Matrix[CSL[r]] ;
         while( p->next != NULL )
         { 
           p = p->next; 
           if(Initial_Flag[p->neighbor]==0) Initial_Matrix[p->neighbor][k]++; 
          }   
         count++;
        }
    }
    
//c. compute the the cost function value of the generated initial solution
   f = 0;
   for(i=0;i<N;i++)
   {
         p = A_Matrix[i];
         while( p->next != NULL )
         { 
           p = p->next; 
           if(patition[p->neighbor] == patition[i]) f++; 
          }   
   }
   
   for(i=0;i<K;i++) SizeG[i] =0;
   for(j=0;j<N;j++) SizeG[patition[j]]++;
   delete [] Inital_Set;

} // To generate an initial solution for k-ECP

/*****************************************************************************/
/*****************  2. One-Swap-Moves Neighborhood  Tabu Search **************/
/*****************************************************************************/ 
//2.1 Clear the delta matrix 
void Clear_Delta_Matrix( )
{
   int x, v ;
   f = 0;
   for( x = 0 ; x < N ; x++ )
     for( v = 0 ; v < K ; v++ )
        Delta_Matrix[ x ][ v ] = 0 ;

   for( x = 0 ; x < N ; x++ )
     for( v = 0 ; v < K ; v++ )
         TabuTenure[ x ][ v ] = 0 ;
   return ;     
}

//2.2 Build delta matrix
void Build_Delta_Matrix( )
{
   int i, j ;
   Clear_Delta_Matrix( ) ;
   for( i = 0 ; i < N ; i++ )
      for( j = 0 ; j < i ; j++ )
         if( Edge[ i ][ j ] != 0 )
           {
             Delta_Matrix[ i ][ Color[ j ] ] += 1;
             Delta_Matrix[ j ][ Color[ i ] ] += 1;
             if( Color[ i ] == Color[ j ] ) 
                 f += 1;
           }
   return ;        
}

//2.3 Update delta matrix after one move
void One_Move_Update_Delta_Matrix(int i, int v0, int v1)
{
   int j ;
   p1 = A_Matrix[ i ] ;
   while( p1->next != NULL )    
        {
          p1 = p1->next ;
          j = p1->neighbor ;
          Delta_Matrix[ j ][ v0 ] -= Edge[ j ][ i ];
          Delta_Matrix[ j ][ v1 ] += Edge[ j ][ i ];
        }                                                                        
   return ;     
}

//2.4 Tabu Search with the union neighborhood of N1 and N2 
int One_swap_Move_Tabu_Search(int Color_into[], int SizeGroup[], int *value, int alpha, int TT)
{
     int tabu_best_delta, best_delta, delt ;
     int num_tabu_best, num_best ;          // the number of tabu neighbors and non-tabu neighbors for the 1-move neighborhood
     int best_x[ 50 ], best_v[ 50 ], x, y, v;
     int tabu_best_x[ 50 ], tabu_best_v[ 50 ];
     
     int swap_num_tabu_best, swap_num_best;  // the number of tabu neighbors and non-tabu neighbors for the swap-move neighborhood
     int swap_best_x[ 50 ],  swap_best_y[ 50 ];
     int swap_tabu_best_x[ 50 ], swap_tabu_best_y[ 50 ];
     int swap_tabu_best_delta, swap_best_delta;
     
     int old_color, old_color1, old_color2;
     int select;
     int num, num2 ;  // records the size of neighborhoods N1 and N2. (the varables (i.e., num, num2) are for only debug and they are not necessary.)
     int cs;   // Records the number of conficting vertices in the current solution. 
     int aver_len, num_conf ;
     int iter;
     int i,j ;
     int non_improve = 0 ;   // The stop condition of TS
     
     int t=0,p=0;
     
     const int p_max = 15; 
     
     T_max = TT ;
          
     int B[p_max] = {1,2,1,4,1,2,1,8,1,2,1,4,1,2,1};
     int A[p_max];
     int T[p_max];
        
     for(i=0;i<p_max;i++)
	 {
         A[i] = int (3.0*T_max*B[i]/8); 
         T[i] = int (1.0*T_max*B[i]/8); 
     }
        
     for(i=0;i<N;i++) Color[i] = Color_into[i]; 
     Build_Delta_Matrix( );
     f_best = f ;
     if(f==0){ *value =0; return 0;}
      
    // printf("\n\n");
    // cout << endl << "One_Move :       iter       f       f_best      time " << endl;
    // cout << "---------------------------------------------------" << endl;
      
     iter = 0 ; 
     while( non_improve < alpha )
        {
          num = 0 ;
          tabu_best_delta = 9999999 ; 
          best_delta = 9999999 ;
          num_tabu_best = 0 ; 
          num_best = 0 ;
          
          num2 = 0; 
          swap_tabu_best_delta = 999999 ; 
          swap_best_delta = 999999 ;
          swap_num_tabu_best = 0 ; 
          swap_num_best = 0 ;
          
          cs=0;
       // a. evaluating the 1-move neighborhood N1
          for( x = 0 ; x < N ; x++ )
            if( Delta_Matrix[ x ][ Color[ x ] ] )
              {
                cs++;
                for( v = 0 ; v < K ; v++ )
                  if( (v != Color[ x ]) && (SizeGroup[Color[x]] > Lbound) && (SizeGroup[v] < Ubound) )
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
                       else if(TabuTenure[ x ][ v ] > iter)// if it is tabu 
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
                 
           // b. evaluating the 2-move neighborhood N2      
            for( x = 0 ; x < N ; x++ )
             if(Delta_Matrix[ x ][Color[x]])
              {
                 for( y = 0 ; y < N ; y++ )
                  if((Color[ x ] !=  Color[ y ]) && (x!=y) )
                    {
						num2 ++ ;
						delt = (Delta_Matrix[ x ][Color[y]] - Delta_Matrix[ x ][ Color[ x ] ]) + (Delta_Matrix[ y ][Color[x]] - Delta_Matrix[ y ][ Color[ y ]]) - 2*Edge[x][y] ;
						if( (TabuTenure[ x ][ Color[y] ] <= iter) && (TabuTenure[ y ][ Color[x] ] <= iter) ) // if this is not tabued 
                        {
							if( delt < swap_best_delta )
							{
								swap_best_x[ 0 ] = x ; 
								swap_best_y[ 0 ] = y ;
								swap_best_delta = delt ; 
								swap_num_best = 1 ;
							}
							else if( delt == swap_best_delta && swap_num_best < 50 )
							{
								swap_best_x[ swap_num_best ] = x ; 
								swap_best_y[ swap_num_best ] = y ;
								swap_num_best ++ ;
							}
						}                                                    
						else if( (TabuTenure[ x ][ Color[y] ] > iter) || (TabuTenure[ y ][ Color[x] ] > iter)) // if it is tabu 
						{ 
							if( delt < swap_tabu_best_delta  )
							{
								swap_tabu_best_x[ 0 ] = x ; 
								swap_tabu_best_y[ 0 ] = y ; 
								swap_tabu_best_delta = delt ; 
								swap_num_tabu_best = 1 ;
							}
							else if( delt == swap_tabu_best_delta && swap_num_tabu_best < 50 )
							{
								swap_tabu_best_x[ swap_num_tabu_best ] = x ; 
								swap_tabu_best_y[ swap_num_tabu_best ] = y ; 
								swap_num_tabu_best ++ ;
							}                               
						}
					}  
                }
           // printf("\n f= %d  |c(s)| = %d   |N1| = %d   |N2| = %d ",f, cs, num, num2); 
           //choose the tabu best move if the tab aspiration criterion is satisfied
           if( (swap_num_tabu_best > 0 && swap_tabu_best_delta < swap_best_delta && ( f + swap_tabu_best_delta < f_best )) || ( num_tabu_best > 0 && tabu_best_delta < best_delta && ( f + tabu_best_delta < f_best ) ) || (num_best + swap_num_best) == 0 )  // aspiration criterion 
                   {
                     
                     if( tabu_best_delta <= swap_tabu_best_delta )
                     {
                       f += tabu_best_delta ;
                       select = rand( ) % num_tabu_best ;  
        
                       old_color = Color[ tabu_best_x[ select ] ] ;              
                       One_Move_Update_Delta_Matrix( tabu_best_x[ select ], old_color, tabu_best_v[ select ] );
                       Color[ tabu_best_x[ select ] ] = tabu_best_v[ select ] ;
					   SizeGroup[old_color]--; SizeGroup[tabu_best_v[ select ]]++;
					   if( iter % 90000 < 30000 )
					   {
                         TabuTenure[ tabu_best_x[ select ] ][ old_color ] = T[p] + rand()%3 ;  
                         TabuTenure[ tabu_best_x[ select ] ][ old_color ] += iter ; 
                       }
                       else if(iter % 90000 > 60000)
                       {
                         TabuTenure[ tabu_best_x[ select ] ][ old_color ] = int (0.9*cs) + rand()%5;  
                         TabuTenure[ tabu_best_x[ select ] ][ old_color ] += iter ;      
                       }
                       else 
                       {
                         TabuTenure[ tabu_best_x[ select ] ][ old_color ] = 5 + rand()%5;  
                         TabuTenure[ tabu_best_x[ select ] ][ old_color ] += iter ;   
                       }
                       
                     }
                     
                     else
                     {
                       f += swap_tabu_best_delta ;
                       select = rand( ) % swap_num_tabu_best ;  
        
                       old_color1 = Color[ swap_tabu_best_x[ select ] ] ;  
                       old_color2 = Color[ swap_tabu_best_y[ select ] ] ; 
                                   
                       One_Move_Update_Delta_Matrix( swap_tabu_best_x[ select ], old_color1, old_color2 );
                       One_Move_Update_Delta_Matrix( swap_tabu_best_y[ select ], old_color2, old_color1 );
                       
                       Color[ swap_tabu_best_x[ select ] ] = old_color2 ;
                       Color[ swap_tabu_best_y[ select ] ] = old_color1 ;
                       
                       if( iter % 90000 < 30000 )
                       {
                         TabuTenure[ swap_tabu_best_x[ select ] ][ old_color1 ] = T[p] + rand()%3 ;  
                         TabuTenure[ swap_tabu_best_y[ select ] ][ old_color2 ] = T[p] + rand()%3 ;
                       }
                       else if(iter % 90000 > 60000)
                       {
                         TabuTenure[ swap_tabu_best_x[ select ] ][ old_color1 ] = int (0.9*cs) + rand()%5;  
                         TabuTenure[ swap_tabu_best_y[ select ] ][ old_color2 ] = int (0.9*cs) + rand()%5; 
                       } 
                       else 
                       {
                         TabuTenure[ swap_tabu_best_x[ select ] ][ old_color1 ] = 5 + rand()%5;  
                         TabuTenure[ swap_tabu_best_y[ select ] ][ old_color2 ] = 5 + rand()%5; 
                       }
                         
                       TabuTenure[ swap_tabu_best_x[ select ] ][ old_color1 ] += iter ; 
                       TabuTenure[ swap_tabu_best_y[ select ] ][ old_color2 ] += iter ;  	
                     }
                     
                     t ++; 
                     
			         if( t > A[p]) {  p = ( p + 1 ) % p_max;   t=0; }	 
                   } 
                   
              else // choose the best non tabu move
                   {
                     if(best_delta <= swap_best_delta)
                     {
                       f += best_delta ; 
                       select = rand( ) % num_best ;    
                       old_color = Color[ best_x[ select ] ] ;              
                       One_Move_Update_Delta_Matrix( best_x[ select ], old_color, best_v[ select ] );
                       Color[ best_x[ select ] ] = best_v[ select ] ;
                       SizeGroup[old_color]--; SizeGroup[best_v[ select ]]++;
                      if( iter % 90000 < 30000 )
                       {
                         TabuTenure[ best_x[ select ] ][ old_color ] = T[p] + rand()%3 ; 
                         TabuTenure[ best_x[ select ] ][ old_color ] += iter ; 
                       }
                       else if(iter % 90000 > 60000)
                       {
                         TabuTenure[ best_x[ select ] ][ old_color ] = int (0.9*cs) + rand()%5; 
                         TabuTenure[ best_x[ select ] ][ old_color ] += iter ;      
                       }
                       else 
                       {
                         TabuTenure[ best_x[ select ] ][ old_color ] = 5 + rand()%5; 
                         TabuTenure[ best_x[ select ] ][ old_color ] += iter ;    
                       }
                       
                     }
                     else 
                     {
                       f += swap_best_delta ;
                       select = rand( ) % swap_num_best ;  
        
                       old_color1 = Color[ swap_best_x[ select ] ] ;  
                       old_color2 = Color[ swap_best_y[ select ] ] ; 
                                   
                       One_Move_Update_Delta_Matrix( swap_best_x[ select ], old_color1, old_color2 );
                       One_Move_Update_Delta_Matrix( swap_best_y[ select ], old_color2, old_color1 );
                       
                       Color[ swap_best_x[ select ] ] = old_color2 ;
                       Color[ swap_best_y[ select ] ] = old_color1 ;
					 
		              if( iter % 90000 < 30000 )
					   {
                       TabuTenure[ swap_best_x[ select ] ][ old_color1 ] = T[p] + rand()%3 ;  
                       TabuTenure[ swap_best_y[ select ] ][ old_color2 ] = T[p] + rand()%3 ; 
                       }
                       else if(iter % 90000 > 60000)
                       {
                         TabuTenure[ swap_best_x[ select ] ][ old_color1 ] = int (0.9*cs) + rand()%5;   
                         TabuTenure[ swap_best_y[ select ] ][ old_color2 ] = int (0.9*cs) + rand()%5;    
                       }
                       else 
                       {
                         TabuTenure[ swap_best_x[ select ] ][ old_color1 ] = 5 + rand()%5;   
                         TabuTenure[ swap_best_y[ select ] ][ old_color2 ] = 5 + rand()%5;     
                       }
                       TabuTenure[ swap_best_x[ select ] ][ old_color1 ] += iter ; 
                       TabuTenure[ swap_best_y[ select ] ][ old_color2 ] += iter ;   
                     } 
                     
                     t++; 
                     
			         if( t > A[p] ) {  p = (p + 1 ) % p_max;    t=0;    }	 
			         
                  } 
                 iter ++ ;
                 //printf("\n One_swap_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
                 total_time = (clock() - starting_time )/CLOCKS_PER_SEC ;
                 if( f <= f_best )
                   {
                     if( f < f_best )
                       {
                         f_best = f;
                         for( i = 0 ; i < N ; i ++ )
                              Best_Color[ i ] = Color[ i ] ;
                        // printf("\n One_swap_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
                         non_improve = 0 ;
                       }  
                     else if ( f == f_best )   
                         non_improve ++ ;
                     if( f_best == 0 )
                        {
                          //printf("\n == %d      %d      %d  ", iter, f, f_best);
                          break; 
                        }  
                   }  
                 else  non_improve ++ ;   
     }
   for(i=0;i<N;i++) Color_into[i] = Best_Color[i];
   for(i=0;i<K;i++) SizeGroup[i]=0; 
   for(i=0;i<N;i++) SizeGroup[Color_into[i]]++; 
   (*value) = f_best; 
   return f_best;
}

/*****************************************************************************/
/*****************          3. Outputing  results      ***********************/
/*****************************************************************************/ 
void Outputing(Solution &S,char *filename)
{
    int i;int r;
	FILE *fp; 
	char buff[80];
	r= rand()%10000;
    sprintf(buff,"%s-results-%d-%d.txt",filename,r, K);
    fp=fopen(buff,"a+");
    fprintf(fp,"N = %d  K = %d  f = %d \n", N , K, S.value); 
    for(i=0;i<K;i++)
    fprintf(fp,"%5.2d \n",S.SizeG[i]); 
    printf("\n");
    for(i=0;i<N;i++)
    fprintf(fp,"%5.2d   %5.2d \n",i, S.p[i]); 
	fclose(fp);
}
void Outresulting(Statistics DATE[], char *filename, char *File_Name, int runs)
{
    int i,j;
    FILE *fp; 
   	char buff[80];
    int N_hit =0;
    double K_SUM =0; 
    double tot_time=0; 
    double time_avg = 0;
    int min_k = 999999;
    sprintf(buff,"%s.txt",filename);
    
    for(i=0; i < runs; i ++)
    { 
      if( DATE[i].K_best <  min_k ) min_k = DATE[i].K_best;
      K_SUM += DATE[i].K_best;
    }
    K_SUM /= runs; 
    
    for(i=0; i < runs; i ++) if( DATE[i].K_best == min_k ) { N_hit ++; tot_time += DATE[i].Time_hit; }
    time_avg = tot_time/N_hit;
    
    fp=fopen(buff,"a+");
    fprintf(fp,"%s   %d    %d    %lf    %d    %lf \n", File_Name, N, min_k, K_SUM, N_hit, time_avg );  
	fclose(fp);         
}
/*****************************************************************************/
/*****************         4. Perturbation operators  *************************/
/*****************************************************************************/ 
void Perturbation(Solution &S, int SizeGroup[], int theta)
{
     int r, delt ;
     int num_best ;      // the number of tabu neighbors and non-tabu neighbors for the 1-move neighborhood
     int best_x[ 50 ], best_v[ 50 ], x, y, v;
     int swap_num_best;  // the number of tabu neighbors and non-tabu neighbors for the swap-move neighborhood
     int swap_best_x[ 50 ], swap_best_y[ 50 ];
     int old_color, old_color1, old_color2;
     int select;
     int num, num2, cs;  // records the size of neighborhoods N1 and N2.
     int iter;
     int i;
   
     for(i=0;i<N;i++) Color[i] = S.p[i]; 
     Build_Delta_Matrix( );
     //printf(" f = %d  \n",f);
     //printf("\n\n");
     //cout << endl << "One_Move :       iter       f       f_best      time " << endl;
     //cout << "-------------------------------------------------------------" << endl;
      
     iter = 0 ; 
     while( iter < theta )
        {
          num2 = 0; 
          swap_num_best = 0;
          cs = 0;                  
       // b. searching the 2-move neighborhood N2  
          for( x = 0 ; x < N ; x++ )
             if(Delta_Matrix[ x ][Color[x]])
              {
                 for( y = 0 ; y < N ; y++ )
                   if((Color[ x ] !=  Color[ y ]) && (x!=y) )
                    {
						num2 ++ ;
					    if( swap_num_best < 50 )
							{
								swap_best_x[ swap_num_best ] = x ; 
								swap_best_y[ swap_num_best ] = y ;
								swap_num_best ++ ;
							}
						else if( swap_num_best == 50)
							{
                                r = rand()% swap_num_best;
								swap_best_x[ r ] = x ; 
								swap_best_y[ r ] = y ;	
							}

					}  
                }
             if(num2==0) break;   
           // printf("\n f= %d  |c(s)| = %d   |N1| = %d   |N2| = %d ",f, cs, num, num2); 
                         
              select = rand( ) % swap_num_best ;  
              x = swap_best_x[ select ];
              y = swap_best_y[ select ]; 
              delt = (Delta_Matrix[ x ][Color[y]] - Delta_Matrix[ x ][ Color[ x ] ]) + (Delta_Matrix[ y ][Color[x]] - Delta_Matrix[ y ][ Color[ y ]]) - 2*Edge[x][y] ;
              f += delt ;
              old_color1 = Color[ x ] ;  
              old_color2 = Color[ y ] ; 
                                   
              One_Move_Update_Delta_Matrix( x, old_color1, old_color2 );
              One_Move_Update_Delta_Matrix( y, old_color2, old_color1 );
                       
              Color[ x ] = old_color2 ;
              Color[ y ] = old_color1 ;   
                                               
              iter ++ ;
             // total_time = (clock() - starting_time )/CLOCKS_PER_SEC ;                                 
     }
   for(i=0;i<N;i++) S.p[i] = Color[i];
   for(i=0;i<K;i++) SizeGroup[i]=0; 
   for(i=0;i<N;i++) SizeGroup[S.p[i]]++; 
   S.value = f; 
}

int Perturb_TS(int Color_into[], int SizeGroup[], int *value, int Iter_max, int L, int *f_cb )
{
     int tabu_best_delta, best_delta, delt ;
     int num_tabu_best, num_best ;  // the number of tabu neighbors and non-tabu neighbors for the 1-move neighborhood
     int best_x[ 50 ], best_v[ 50 ], x, y, v;
     int tabu_best_x[ 50 ], tabu_best_v[ 50 ];
     
     int swap_num_tabu_best, swap_num_best;  // the number of tabu neighbors and non-tabu neighbors for the swap-move neighborhood
     int swap_best_x[ 50 ],  swap_best_y[ 50 ];
     int swap_tabu_best_x[ 50 ], swap_tabu_best_y[ 50 ];
     int swap_tabu_best_delta, swap_best_delta;
     
     int old_color, old_color1, old_color2;
     int select;
     int num, num2 ; // records the size of neighborhoods N1 and N2.
     int cs;  // Records the size of set of conficting vertices. (the varables (i.e., num, num2, and cs) are for only debug and they are not necessary.)
     int aver_len, num_conf ;
     int iter;
     int i,j ;
     int non_improve = 0 ;  // The stop condition of TS
     
     int t=0,p=0;
     const int p_max = 15; 
     int TT_max = L ; 
         
     for(i=0;i<N;i++) Color[i] = Color_into[i]; 
     for(i=0;i<N;i++) Best_Color[i] = Color_into[i];
     Build_Delta_Matrix();
     f_best = f ;
     //printf("\n\n");
     //cout << endl << "One_Move :       iter       f       f_best      time " << endl;
     //cout << "-------------------------------------------------------------" << endl;
      
     iter = 0 ; 
     while( iter < Iter_max )
        {
          num = 0 ;
          tabu_best_delta = 9999999 ; 
          best_delta = 9999999 ;
          num_tabu_best = 0 ; 
          num_best = 0 ;
          
          num2 = 0; 
          swap_tabu_best_delta = 999999 ; 
          swap_best_delta = 999999 ;
          swap_num_tabu_best = 0 ; 
          swap_num_best = 0 ;
          
          cs=0;
       // a. evaluating the 1-move neighborhood N1
          for( x = 0 ; x < N ; x++ )
            if( Delta_Matrix[ x ][ Color[ x ] ] )
              {
                cs++;
                for( v = 0 ; v < K ; v++ )
                  if( (v != Color[ x ]) && (SizeGroup[Color[x]] > Lbound) && (SizeGroup[v] < Ubound) )
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
                       else if(TabuTenure[ x ][ v ] > iter)// if it is tabu 
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
                 
           // b. evaluating the 2-move neighborhood N2      
            for( x = 0 ; x < N ; x++ )
             if(Delta_Matrix[ x ][Color[x]])
              {
                 for( y = 0 ; y < N ; y++ )
                  if((Color[ x ] !=  Color[ y ]) && (x!=y) )
                    {
						num2 ++ ;
						delt = (Delta_Matrix[ x ][Color[y]] - Delta_Matrix[ x ][ Color[ x ] ]) + (Delta_Matrix[ y ][Color[x]] - Delta_Matrix[ y ][ Color[ y ]]) - 2*Edge[x][y] ;
						if( (TabuTenure[ x ][ Color[y] ] <= iter) && (TabuTenure[ y ][ Color[x] ] <= iter) ) // if this is not tabued 
                        {
							if( delt < swap_best_delta )
							{
								swap_best_x[ 0 ] = x ; 
								swap_best_y[ 0 ] = y ;
								swap_best_delta = delt ; 
								swap_num_best = 1 ;
							}
							else if( delt == swap_best_delta && swap_num_best < 50 )
							{
								swap_best_x[ swap_num_best ] = x ; 
								swap_best_y[ swap_num_best ] = y ;
								swap_num_best ++ ;
							}
						}                                                    
						else if( (TabuTenure[ x ][ Color[y] ] > iter) || (TabuTenure[ y ][ Color[x] ] > iter)) // if it is tabu 
						{ 
							if( delt < swap_tabu_best_delta  )
							{
								swap_tabu_best_x[ 0 ] = x ; 
								swap_tabu_best_y[ 0 ] = y ; 
								swap_tabu_best_delta = delt ; 
								swap_num_tabu_best = 1 ;
							}
							else if( delt == swap_tabu_best_delta && swap_num_tabu_best < 50 )
							{
								swap_tabu_best_x[ swap_num_tabu_best ] = x ; 
								swap_tabu_best_y[ swap_num_tabu_best ] = y ; 
								swap_num_tabu_best ++ ;
							}                               
						}
					}  
                }
           // printf("\n f= %d  |c(s)| = %d   |N1| = %d   |N2| = %d ",f, cs, num, num2); 
           //choose the tabu best move if the tab aspiration criterion is satisfied
           if( (swap_num_tabu_best > 0 && swap_tabu_best_delta < swap_best_delta && ( f + swap_tabu_best_delta < f_best )) || ( num_tabu_best > 0 && tabu_best_delta < best_delta && ( f + tabu_best_delta < f_best ) ) || (num_best + swap_num_best) == 0 )  // aspiration criterion 
                   {
                     
                     if( tabu_best_delta <= swap_tabu_best_delta )
                     {
                       f += tabu_best_delta ;
                       select = rand( ) % num_tabu_best ;  
        
                       old_color = Color[ tabu_best_x[ select ] ] ;              
                       One_Move_Update_Delta_Matrix( tabu_best_x[ select ], old_color, tabu_best_v[ select ] );
                       Color[ tabu_best_x[ select ] ] = tabu_best_v[ select ] ;
					   SizeGroup[old_color]--; SizeGroup[tabu_best_v[ select ]]++;
					  
                       TabuTenure[ tabu_best_x[ select ] ][ old_color ] = TT_max + rand()%1000 ;   
                       TabuTenure[ tabu_best_x[ select ] ][ old_color ] += iter ; 

                     }
                     
                     else
                     {
                       f += swap_tabu_best_delta ;
                       select = rand( ) % swap_num_tabu_best ;  
        
                       old_color1 = Color[ swap_tabu_best_x[ select ] ] ;  
                       old_color2 = Color[ swap_tabu_best_y[ select ] ] ; 
                                   
                       One_Move_Update_Delta_Matrix( swap_tabu_best_x[ select ], old_color1, old_color2 );
                       One_Move_Update_Delta_Matrix( swap_tabu_best_y[ select ], old_color2, old_color1 );
                       
                       Color[ swap_tabu_best_x[ select ] ] = old_color2 ;
                       Color[ swap_tabu_best_y[ select ] ] = old_color1 ;
        
                       TabuTenure[ swap_tabu_best_x[ select ] ][ old_color1 ] = TT_max + rand()%1000 ;   
                       TabuTenure[ swap_tabu_best_y[ select ] ][ old_color2 ] = TT_max + rand()%1000 ; 
                       
                       TabuTenure[ swap_tabu_best_x[ select ] ][ old_color1 ] += iter ; 
                       TabuTenure[ swap_tabu_best_y[ select ] ][ old_color2 ] += iter ;  	
                     }
                        
                   } 
              else // choose the best non tabu move
                   {
                     if(best_delta <= swap_best_delta)
                     {
                       f += best_delta ; 
                       select = rand( ) % num_best ;    
                       old_color = Color[ best_x[ select ] ] ;              
                       One_Move_Update_Delta_Matrix( best_x[ select ], old_color, best_v[ select ] );
                       Color[ best_x[ select ] ] = best_v[ select ] ;
                       SizeGroup[old_color]--; SizeGroup[best_v[ select ]]++;
                     
                       TabuTenure[ best_x[ select ] ][ old_color ] = TT_max + rand()%1000 ;  
                       TabuTenure[ best_x[ select ] ][ old_color ] += iter ; 

                     }
                     else 
                     {
                       f += swap_best_delta ;
                       select = rand( ) % swap_num_best ;  
        
                       old_color1 = Color[ swap_best_x[ select ] ] ;  
                       old_color2 = Color[ swap_best_y[ select ] ] ; 
                                   
                       One_Move_Update_Delta_Matrix( swap_best_x[ select ], old_color1, old_color2 );
                       One_Move_Update_Delta_Matrix( swap_best_y[ select ], old_color2, old_color1 );
                       
                       Color[ swap_best_x[ select ] ] = old_color2 ;
                       Color[ swap_best_y[ select ] ] = old_color1 ;
					 
		             
                       TabuTenure[ swap_best_x[ select ] ][ old_color1 ] = TT_max + rand()%1000 ;  
                       TabuTenure[ swap_best_y[ select ] ][ old_color2 ] = TT_max + rand()%1000 ;  
                       
                       TabuTenure[ swap_best_x[ select ] ][ old_color1 ] += iter ; 
                       TabuTenure[ swap_best_y[ select ] ][ old_color2 ] += iter ;   
                     } 
                        
                  } 
                 iter ++ ;
                
                 //printf("\n per_swap_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
                 total_time = (clock() - starting_time )/CLOCKS_PER_SEC ;
                 
                 if( f <= f_best )
                   {
                     if( f < f_best )
                       {
                         f_best = f;
                         for( i = 0 ; i < N ; i ++ )
                         Best_Color[ i ] = Color[ i ] ;
                         // printf("\n Perturb_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );
                         non_improve = 0 ;
                       }  
                     else if ( f == f_best )   
                         non_improve ++ ;
                     if( f_best == 0 )
                        {
                         // printf("\n == %d      %d      %d  ", iter, f, f_best);
                          break; 
                        }  
                   }  
                 
                 else  non_improve ++ ;   
                //  printf("\n Perturb_Move :  %8d       %3d       %3d       %5.3lf s", iter, f, f_best, total_time );   
     }
   (*f_cb) = f_best; 
   for(i=0;i<N;i++) Color_into[i] = Color[i];
   for(i=0;i<K;i++) SizeGroup[i]=0; 
   for(i=0;i<N;i++) SizeGroup[Color[i]] ++; 
   (*value) = f; 
   return f;  
}
/*****************************************************************************/
/**************** 5. Binary Search for Initial value of K ********************/
/*****************************************************************************/ 
void BinarySearch(Solution &S, int *VK)
{    
     int i; 
     int LowerK = 0;
     int UpperK = N;
     int Mid = N/2; 
     G_K = N ;
     while(UpperK > LowerK + 1)
     {
        K = Mid;
        NewBound();             
        Initia_sol2(S.p, S.SizeG); S.value = 99999;
        One_swap_Move_Tabu_Search(S.p, S.SizeG, &S.value, 100, 150); 
        //printf("\n K=%d \n",K);
        if(S.value == 0)
        {
         UpperK =  K ;
         G_K = K; 
         Time_one_run_hit = (clock() - Time_one_run_stating) / CLOCKS_PER_SEC;
         Mid = (LowerK + UpperK)/2 ; 
         for(i=0;i<N;i++) GS_best.p[i] = S.p[i];
         GS_best.k1 = K; 
         //Outputing(S,File_Name);
        }
        else 
        {
          LowerK = Mid ;
          Mid = (LowerK + UpperK)/2 ;
        }
        // free_memery(); 
        //printf("\n lowerK= %d  upperK = %d \n",LowerK,UpperK);
     }  
     
    *VK = UpperK;  
       // printf("K=%d \n",*VK);
}

/*****************************************************************************/
/************* 6. Backtrack based Iterated Tabu Search ***********************/
/*****************************************************************************/ 
void ITS(Solution &S)
{    
     int i,j;
     int f_cb ;
     int f_gb = 99999; 
     int TT = 80, TTmax = 2000;  //It may be better to set TT to 50. 

     int L0 = 5000, LL;
     int Delta_L = 500; 
     double Rand;
     int Npert; 
     int NoImprove; 
     

     Npert = 30;
     
     LL= L0;
     NoImprove = 0;  
     Initia_sol2(S.p, S.SizeG); S.value = 999999; 
     One_swap_Move_Tabu_Search(S.p, S.SizeG, &S.value, ALPHA, TT); 
     if(S.value==0)
      { 
        G_K = K; 
        for(i=0;i<N;i++) GS_best.p[i] = S.p[i]; GS_best.k1 = K; 
        Time_one_run_hit = (clock() - Time_one_run_stating) / CLOCKS_PER_SEC;
        return;
      }
      
     while( (NoImprove < Npert) && ( (clock() - Time_one_run_stating) / CLOCKS_PER_SEC < time_limit ) )
     {
      //printf("\n K=%d \n",K);
      
      for(i=0;i<N;i++) SC.p[i] = S.p[i];
      for(j=0;j<K;j++) SC.SizeG[j] = S.SizeG[j];
      SC.value = S.value; 
      
      Rand = 1.0*(rand()%30000)/30000.0;
      if(Rand < 0.7) Perturbation(SC, SC.SizeG, int(0.3*N) ); 
      else Perturb_TS(SC.p, SC.SizeG, &SC.value, L0, TTmax, &f_cb);      
         
      if(SC.value==0)
       { 
        G_K = K; 
        for(i=0;i<N;i++) GS_best.p[i] = SC.p[i]; GS_best.k1 = K; 
        Time_one_run_hit = (clock() - Time_one_run_stating) / CLOCKS_PER_SEC;
        break;
       }      
       
      One_swap_Move_Tabu_Search(SC.p, SC.SizeG, &SC.value, ALPHA, TT); 
    
      if(SC.value==0)
      { 
        G_K = K; 
        for(i=0;i<N;i++) GS_best.p[i] = SC.p[i]; GS_best.k1 = K; 
        Time_one_run_hit = (clock() - Time_one_run_stating) / CLOCKS_PER_SEC;
        break;
      }
      
      if(SC.value <= S.value)
      {   
          if(SC.value < S.value)
          {
            for(i=0;i<N;i++) S.p[i] = SC.p[i];
            for(j=0;j<K;j++) S.SizeG[j] = SC.SizeG[j];
            S.value = SC.value; 
            NoImprove = 0;
          }
          else NoImprove++; 
           
      }
      else NoImprove++; 
      
    }  
     
}

void Backtrack_ITS( Solution &S )
{   
    int K_0;  
    BinarySearch(S,&K_0);
    K = K_0;
    //printf("K = %d ********************************** \n",K);
    while(1)
     {    
         if( K <= G_K - 1 || K <= 3) K = G_K - 1; 
         else K--;
         NewBound();
         ITS(S);
         if( (clock() - Time_one_run_stating) / CLOCKS_PER_SEC > time_limit ) break;
     }
}

int proof(Solution &GS)
{
    int i,j; int K1;
    int flag = 0,f = 0; 
    int L,U; 
    int *Size;
    K1 = GS.k1 ;
    Size = new int [K1]; 
    L =  int (floor(1.0*N/K1));
    U =  int (ceil(1.0*N/K1)) ; 
    
    for(i=0; i<N; i++)
      for(j=i+1; j<N; j++)
       if(Edge[i][j]==1 && GS.p[i] == GS.p[j]) f++; 
     
    for(i=0;i<K1;i++) Size[i] = 0 ;
    for(i=0;i<N;i++) Size[GS.p[i]]++; 
    
    for(i=0;i<K1;i++)  
      if(Size[i] < L || Size[i] > U ) {flag = 1; break;}
      
    if( flag == 0 && f == 0) 
    {   
        printf("\n K = %d  f = %d \n",K1,f);
        for(i=0;i<K1;i++) printf("%d ",Size[i]); 
        printf("\n");
        for(j=0;j<N;j++) printf("%d ",GS.p[j]); 
    }
    delete [] Size; Size= NULL; 
}
/*****************************************************************************/
/*****************          7. Main Scheme            ************************/
/*****************************************************************************/ 
int main(int argc, char **argv)
{ 
     int i, j, seed ; 
     int sum=0;
     int N_runs = 20;
     int x;
     starting_time = clock(); 
     seed = time(NULL) % 100000 ;
     srand( seed ) ;
     /*
     if( argc == 4 )
       { 
         File_Name = argv[1];
         K = atoi(argv[2]);
         time_limit = atof(argv[3]);
         seed = time(NULL) % 1000000 ;
         srand( seed ) ;
       }
     else 
       {
           cout << endl << "error :" << endl;
           cout << "ITS file_name K time_limit" << endl;
           exit(0);
       }
     */
     // File_Name = "qg.order60.col";
     // outfilename = "ss.txt";
    
     File_Name = argv[1];
     outfilename = argv[2];

     outfilename = "ss1.txt";

     Initializing();
     if(N<=500) time_limit = 10000;
     else time_limit = 20000;
     S.p = new int [N];
     S.SizeG = new int [N/2];
     SC.p = new int [N];
     SC.SizeG = new int [N/2];
     GS_best.p = new int [N];
     GS_best.SizeG = new int [N/2];    
     for(x=0;x<N;x++) Delta_Matrix[x]= new int[N/2];
     for(x=0;x<N;x++) Initial_Matrix[x]= new int[N/2];
     for(x=0;x<N;x++) TabuTenure[x]= new int[N/2];

     for(i=0;i<N_runs;i++)
     { 
       G_K = 999999; 
       Time_one_run_stating = clock();
       Backtrack_ITS(S);
       Results[i].Time_hit = Time_one_run_hit;
       Results[i].K_best = G_K;
       proof(GS_best);
     }
     
     for(i=0;i<N_runs;i++) printf("\n K_BEST = %d  Time = %lf \n", Results[i].K_best, Results[i].Time_hit);
     Outresulting(Results, outfilename, File_Name, N_runs);
     total_time = (clock() - starting_time )/CLOCKS_PER_SEC ; 
     printf("Total time = %lf Seconds \n",total_time); 
     free_memery();
     return 0;
}
 
