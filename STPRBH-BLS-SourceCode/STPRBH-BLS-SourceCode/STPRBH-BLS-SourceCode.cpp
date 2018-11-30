#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <string.h>
#include <math.h>

using namespace std;
#define Null -1
#define InfiniteEdgeCost  1000000
#define MaxInstNum 400
#define MaxVtxNum 500
#define MaxPathLength    50
#define MaxStoredSolutionNum 10000
#define RunBLSTimes 10

char *File_Name;
char *InputInstFile = "A-InstancesToRun-C2.txt";
int Series=3;   //1<-->B  2<-->C1 3<-->C2 4<-->C3 5<-->C4
char OutputRltFile[100];
char inputfilename[100];
char outfilename[100];
int InputParam[MaxInstNum][5];

int InstNum;
int *VtxRevenue;
int *VtxWithPlusRevenue;
int *VtxOnSolution;
int *SimplyVtxIndex;
bool *IfRemoveVtx;
int **EdgeBudget;      
int **TempEdgeBudget;
int **AdjacentVtx;
int **ShortestPath;
int ***CompleteShortestPath;
int **TempShortestPath;

int **ShortestPathMatrix;
int *ShortestPathToCurSlt;
int *TempShortestPathToCurSlt;

int *FinalSolution;
int *TempSolutionArray;
int *CurIterBestSolutionArray;
int *BestSolutionArray;
int *BLSBestSolutionArray;
int *GlobalBestSolutionArray;
int *BackBone;
int *RetainedSolution;

int **RecordSolution;
int **StoredSolution;
int *StoredSolutionHashValue;
int StoredIndex;
int *VtxOccurTimes;

bool *IfVtxOnSolution;
int  *VtxStatus;
bool *IfLeafVtx;
float *VtxPriority;

int Max_Vtx; 
int Max_Edge;
int Max_Budget;
int VtxWithPlusRvnNum;
int ExpectedRevenue;
int SltVtxNum;
int SltEdgNum;
int CurSltBudget;
int ProfitNodeNum;
int RecordPath[MaxPathLength];
int PathLength;
int MaxHop=13;
int BudgetDivideNum=10;

int CandidateVtx[MaxVtxNum];
float VtxFitniss[MaxVtxNum];
int ShortestPathConnectedVtx[MaxVtxNum];
int VtxWithFitNum=0;

// Parameters for running the BLS
int Version=0;
int Alpha=Null;
int MaxRecordSolutionNum=100; 
int PrabToAddBestVtx=30;
int MaxNS=3;
double LimitedDis=0.3;
int MaxNoImpPerturb=100;
double MaxUsedTime=720000;

int DisplayMode=Null; 
bool IfDisplayDetail;
int RecordSltNum;
int AvgDistance;
double starting_time;
double finishing_time;
double avg_time;

int InstCls;
int InstNo;
int InstType;
double PrepareTime;
int CollectedRev[RunBLSTimes];
double UseTime[RunBLSTimes];

int RandomInt(int DivideNum)
{
   int RandomNum =rand() % DivideNum;
   return RandomNum;
}// End RandomInt

int GetRootNode()
{
   for(int i=0; i<Max_Vtx; i++)
   {
     if(VtxRevenue[i]>0)
       return i;        
   }    
}//End GetRootNode()

void GetAdjacentMatrix()
{
   for(int i=0;i<Max_Vtx;i++)
     for (int j=0;j<Max_Vtx;j++)        
       AdjacentVtx[i][j]=0;
       
   for(int i=0;i<Max_Vtx-1;i++)
     for(int j=i+1;j<Max_Vtx;j++) 
     {
        if(EdgeBudget[i][j] != InfiniteEdgeCost)
        {
           AdjacentVtx[i][0]++;
           AdjacentVtx[i][AdjacentVtx[i][0]]=j;
           AdjacentVtx[j][0]++;
           AdjacentVtx[j][AdjacentVtx[j][0]]=i;                              
        }//end if
     }//end for     
}//end GetAdjacentMatrix()

void SwapVertex(int FirstVtx, int SecondVtx)
{
   if(FirstVtx==SecondVtx)
     return;
     
   int TempRevenue;
   TempRevenue =  VtxRevenue[FirstVtx];
   VtxRevenue[FirstVtx] = VtxRevenue[SecondVtx];
   VtxRevenue[SecondVtx] = TempRevenue;
   
   int TempCost;
   for(int i=0;i<Max_Vtx;i++)
   {
      if(i!=FirstVtx && i!=SecondVtx)
      {
          TempCost = EdgeBudget[FirstVtx][i];
          EdgeBudget[FirstVtx][i] = EdgeBudget[SecondVtx][i];
          EdgeBudget[SecondVtx][i] = TempCost;
     
          TempCost = EdgeBudget[i][FirstVtx];
          EdgeBudget[i][FirstVtx] =EdgeBudget[i][SecondVtx];
          EdgeBudget[i][SecondVtx] = TempCost;
       }
   }
   
   for(int i=0;i<Max_Vtx;i++)
     for(int j=0; j< Max_Vtx; j++)
     TempEdgeBudget[i][j]=EdgeBudget[i][j];
        
   GetAdjacentMatrix();
}//End SwapVertex()

//Get the total revenue of the original graph
int GetTotalRevenue()
{
   int Total_Revenue=0;          
   for(int i=0;i<Max_Vtx;i++)  
      Total_Revenue += VtxRevenue[i];                 

   return Total_Revenue;    
}//End GetTotalRevenue

//Get the total revenue of the current solution
int GetSolutionRevenue()
{    
   int Solution_Revenue=0;      

    for(int i=0; i<VtxWithPlusRvnNum; i++) 
    {              
       int CurVtx=VtxWithPlusRevenue[i];
       if(IfVtxOnSolution[CurVtx]) 
         Solution_Revenue += VtxRevenue[CurVtx];  
    }             

   return Solution_Revenue;    
}//End GetTotalRevenue

//Get the total budget of the original graph
int GetTotalBudget()
{
   int Total_Budget=0;          
   for(int i=0;i<Max_Vtx-1;i++)
     for(int j=i+1;j<Max_Vtx;j++)
     {
        if(TempEdgeBudget[i][j]<InfiniteEdgeCost)
           Total_Budget += TempEdgeBudget[i][j];                 
     }
    //cout<<endl<<"The total budget of the current solution is: "<<Total_Budget<<endl;
    return Total_Budget;    
}//End GetTotalBudget

//Get the total budget of the current solution
int GetSolutionBudget()
{
   int Solution_Budget=0;           //The total budget of the current solution
   for(int i=0; i<Max_Vtx;i++)
   {
      if(FinalSolution[i]!=Null) 
        Solution_Budget += TempEdgeBudget[FinalSolution[i]][i];       
   }   
   //cout<<endl<<"The total budget of the current solution is: "<<Total_Budget<<endl;
   return Solution_Budget;    
}//End GetSolutionBudget

void FloydWarshall()
{
   for(int i=0; i<Max_Vtx; i++)
     for(int j=0; j<Max_Vtx; j++)    
       ShortestPathMatrix[i][j]=EdgeBudget[i][j];  
  
   for(int k=0; k<Max_Vtx; k++)
     for(int i=0; i<Max_Vtx; i++)
       for(int j=0; j<Max_Vtx; j++) 
         if(ShortestPathMatrix[i][k] + ShortestPathMatrix[k][j] < ShortestPathMatrix[i][j])
            ShortestPathMatrix[i][j] = ShortestPathMatrix[i][k] + ShortestPathMatrix[k][j]; 
}//End FloydWarshall()

void Initializing()
{
   srand( (unsigned)time(NULL) ); 
        
   ifstream FIC;
   FILE *fp;   
   FIC.open(inputfilename);   
   if (FIC.fail())
   {
     cout << "### Fail to open, File_Name " << File_Name << endl;
     getchar();
     exit(0);
   }
   char StrReading[100];     
   FIC >> StrReading;
   if (FIC.eof() )
   {
     cout << "### Wrong begin symbol, File_Name " << File_Name << endl;
     exit(0);
   }
  
   //Initilize
   int TempData;      
   FIC>>Max_Vtx>>Max_Edge>>ProfitNodeNum>>TempData;  
   if(IfDisplayDetail)
   {
     cout << "Initializing "<<endl;    
     cout << "Number of vectices = " <<Max_Vtx<<endl;
     cout << "Number of edges = " <<Max_Edge<<endl;
     cout << "Number of profit nodes = " <<ProfitNodeNum<<endl;  
   }   
          
   VtxRevenue = new int [Max_Vtx];
   VtxWithPlusRevenue = new int [Max_Vtx];
   VtxOnSolution = new int [Max_Vtx];
   SimplyVtxIndex = new int [Max_Vtx]; 
   IfRemoveVtx = new bool [Max_Vtx];                
   EdgeBudget = new int*[Max_Vtx];
   TempEdgeBudget = new int*[Max_Vtx];
   AdjacentVtx = new int*[Max_Vtx];
   ShortestPath = new int*[Max_Vtx];
   CompleteShortestPath = new int **[Max_Vtx];
   ShortestPathToCurSlt = new int [Max_Vtx];
   TempShortestPathToCurSlt = new int [Max_Vtx];
   
   TempShortestPath = new int*[Max_Vtx];
   ShortestPathMatrix = new int*[Max_Vtx];
   VtxOccurTimes = new int[Max_Vtx];      
     
   FinalSolution = new int[Max_Vtx];
   TempSolutionArray = new int[Max_Vtx];
   CurIterBestSolutionArray = new int[Max_Vtx];
   BestSolutionArray = new int[Max_Vtx];
   BLSBestSolutionArray = new int[Max_Vtx];
   GlobalBestSolutionArray = new int[Max_Vtx];
   BackBone = new int[Max_Vtx]; 
   RetainedSolution = new int[Max_Vtx];  
   
   RecordSolution = new int *[MaxRecordSolutionNum];
   for(int i=0; i<MaxRecordSolutionNum; i++)
     RecordSolution[i] = new int [Max_Vtx];
         
   IfVtxOnSolution = new bool [Max_Vtx];
   VtxStatus       = new int  [Max_Vtx];
   IfLeafVtx       = new bool [Max_Vtx];
   VtxPriority     = new float [Max_Vtx];
   
   for( int i=0; i <Max_Vtx ; i++ ) 
     CompleteShortestPath[i]= new int*[Max_Vtx];
     
   for( int i=0; i <Max_Vtx ; i++ )
     for( int j=0; j<Max_Vtx; j++ )
       CompleteShortestPath[i][j]=new int [MaxPathLength];     
          
   for( int i=0; i <Max_Vtx ; i++ ) 
   {
      VtxRevenue[i] = 0; 
      VtxOnSolution[i] = Null;
      SimplyVtxIndex[i]=i;
      IfRemoveVtx[i]=false;        
      EdgeBudget[i] = new int[Max_Vtx];        
      TempEdgeBudget [i] = new int[Max_Vtx];
      AdjacentVtx[i] = new int[Max_Vtx];        
      ShortestPath[i]=new int[MaxPathLength];
      TempShortestPath[i]=new int [MaxPathLength];
      ShortestPathMatrix[i] = new int [Max_Vtx];
      ShortestPathToCurSlt[i] = Null;
      TempShortestPathToCurSlt[i] = Null;
        
      FinalSolution[i]=Null; 
      TempSolutionArray[i]=Null;
      CurIterBestSolutionArray[i]=Null;
      BestSolutionArray[i]=Null;
      BLSBestSolutionArray[i]=Null;
      GlobalBestSolutionArray[i]=Null; 
      BackBone[i]=Null;              
        
      IfVtxOnSolution[i] = false;
      VtxStatus[i]=0;
      IfLeafVtx[i] = false;    
      VtxPriority[i] = 0; 
      VtxOccurTimes[i] = 0;  
   } 
   
   StoredSolution = new int *[MaxStoredSolutionNum];
   StoredSolutionHashValue = new int [MaxStoredSolutionNum];
   StoredIndex=0;
   for(int i=0;i<MaxStoredSolutionNum;i++)
   {    
     StoredSolution[i]=new int [Max_Vtx]; 
     StoredSolutionHashValue[i]=0;
   }
     
   for(int i=0;i<MaxStoredSolutionNum;i++)
     for(int j=0;j<Max_Vtx;j++)
     StoredSolution[i][j]=Null; 
            
   for (int i=0; i<Max_Vtx; i++ )
     for (int j=0; j<Max_Vtx; j++ )
     {          
        if(i==j)
           EdgeBudget[i][j] = 0;
        else
           EdgeBudget[i][j] = InfiniteEdgeCost;       
     }   
     
   //Read the revenue of each vetex       
   int *ProfitNode;
   int TempProfitNode;
   ProfitNode = new int [2*ProfitNodeNum];  
   for(int i=0;i<2*ProfitNodeNum;i++)
   {
      FIC>>TempProfitNode;
      ProfitNode[i]=TempProfitNode;    
   }
     
   for(int i=0;i<ProfitNodeNum;i++)
   {          
      VtxRevenue[ProfitNode[i]-1]=ProfitNode[i+ProfitNodeNum];  
   }   

   int x1,x2,EdgeCost;            
   for(int i=0;i<Max_Edge;i++)
   {
      FIC>>x1>>x2>>EdgeCost;  
      x1--; x2--;
      if ( x1<0 || x2<0 || x1>=Max_Vtx || x2 >=Max_Vtx )
      {
         cout << "### Error of node : x1="<< x1 << ", x2=" << x2 << endl;
         getchar();
         exit(0);
      }
      else
        EdgeBudget[x1][x2]=EdgeBudget[x2][x1]=EdgeCost*100;        
   }
  
   int IfInputSimplyIndex=0;  
   if(IfInputSimplyIndex==1)
   {
      for( int i=0; i<Max_Vtx; i++ ) 
      {
         FIC>>x1>>x2;
         SimplyVtxIndex[x1-1]=x2-1;
       }
    }
    else
    {
       for( int i=0; i <Max_Vtx ; i++ )       
          SimplyVtxIndex[i]=i;          
    }

   if(GetRootNode()!=0)
   {
      cout << "Swap node 1 and node " << GetRootNode()+1 <<endl;
      //Swap node 0 with the first node with revenue>0 
      SwapVertex(0, GetRootNode());
   }
   
   VtxWithPlusRvnNum=0;  
   for(int i=0; i<Max_Vtx; i++)
   {
      if(VtxRevenue[i]>0)
        VtxWithPlusRevenue[VtxWithPlusRvnNum++]=i;     
   }
      
   //Store the cost of each edge to a temp matrix  
   for (int i=0; i<Max_Vtx; i++ )
     for (int j=0; j<Max_Vtx; j++ )
       TempEdgeBudget[i][j] = EdgeBudget[i][j];
   
   //Get the adjacent vetices of each vetex, stored in the matirx AdjacentVtx[][]   
   GetAdjacentMatrix();
   FloydWarshall();       
   FIC.close();  
   if(IfDisplayDetail)
   {
     cout<<"\nThe total revenue of the graph is:"<< GetTotalRevenue()<<endl;
     cout<<"The total cost of the graph is:"<< GetTotalBudget()/100<<endl;     
     cout<<"Initializing finished, success! Input any key to continue:"<<endl;
   }
} // End Initializing

void ReleaseMemory()
{     
   delete []VtxRevenue;   
   delete []VtxWithPlusRevenue;
   delete []VtxOnSolution;
   delete []SimplyVtxIndex;
   delete []IfRemoveVtx;
   
   for(int i=0; i<Max_Vtx;i++)
     delete []EdgeBudget[i];
   delete []EdgeBudget;
   
   for(int i=0; i<Max_Vtx;i++)
     delete []TempEdgeBudget[i];
   delete []TempEdgeBudget;
   
   for(int i=0; i<Max_Vtx;i++)
     delete []AdjacentVtx[i];
   delete []AdjacentVtx;
   
   for(int i=0; i<Max_Vtx;i++)
     delete []ShortestPath[i];
   delete []ShortestPath;
   
   for(int i=0; i<Max_Vtx;i++)
   {
     for(int j=0; j<Max_Vtx;j++)
       delete []CompleteShortestPath[i][j];
     delete []CompleteShortestPath[i];
   }     
   delete []CompleteShortestPath;
  
   delete []ShortestPathToCurSlt;
   delete []TempShortestPathToCurSlt;
   
   for(int i=0; i<Max_Vtx;i++)
     delete []TempShortestPath[i];
   delete []TempShortestPath;
  
   for(int i=0; i<Max_Vtx;i++)
     delete []ShortestPathMatrix[i];
   delete []ShortestPathMatrix;   

   delete []VtxOccurTimes;
   delete []FinalSolution;
   delete []TempSolutionArray;
   delete []CurIterBestSolutionArray;
   delete []BestSolutionArray;
   delete []BLSBestSolutionArray;
   delete []GlobalBestSolutionArray;
   delete []BackBone;
   delete []RetainedSolution;
  
   for(int i=0; i<MaxRecordSolutionNum;i++)
     delete []RecordSolution[i];
   delete []RecordSolution;
  
   for(int i=0; i<MaxStoredSolutionNum;i++)
     delete []StoredSolution[i];
   delete []StoredSolution;   
   
   delete []IfVtxOnSolution;
   delete []VtxStatus;
   delete []IfLeafVtx;
   delete []VtxPriority;

   cout<<"Finished releasing memory"<<endl<<endl<<endl<<endl;
} //End ReleaseMemory()

//Update that if some vertex is on the solution or not
void UpdateSolution()
{  
   IfVtxOnSolution[0]=true; 
   IfLeafVtx[0]=true;
   VtxOnSolution[0]=1;  
   VtxOnSolution[1]=0; 
   SltEdgNum=0; 
   SltVtxNum=1;
   CurSltBudget=0;
           
   for(int i=1; i<Max_Vtx;i++)
   {
      IfVtxOnSolution[i]=false;  
      IfLeafVtx[i]=false;
      
      if(FinalSolution[i]!=Null)
      {
         IfVtxOnSolution[i]=true;          
         IfLeafVtx[i]=true;
         SltEdgNum++;  
         SltVtxNum++; 
         VtxOnSolution[0]++;        
         VtxOnSolution[VtxOnSolution[0]]=i;  
         CurSltBudget += TempEdgeBudget[FinalSolution[i]][i];                                  
      }      
   }
   
   for(int i=2;i<=VtxOnSolution[0];i++)
   {     
     IfVtxOnSolution[FinalSolution[VtxOnSolution[i]]]=true;
     IfLeafVtx[FinalSolution[VtxOnSolution[i]]]=false;
   } 
}//End UpdateSolution

int GetVtxDepth(int VtxNum)
{
   int VtxDepth=0;
   int CurrentVtx=VtxNum;
   int ParentVtx=FinalSolution[CurrentVtx];
   while(ParentVtx!=Null)
   {
      VtxDepth++;
      CurrentVtx=ParentVtx; 
      ParentVtx=FinalSolution[CurrentVtx];                   
   }   
   
   if(CurrentVtx==0)
     return VtxDepth;
   else
     return Null;
}//End GetVtxDepth

//Check if a solution if feasible or unfeasible
//The solution is recorded as a directed tree
bool CheckSolutionFeasible()
{ 
   UpdateSolution();
     
   if(IfDisplayDetail) 
     cout<<"\n Begin to check if the current solution is feasible"<<endl;
   
   //Check if the number of edges is 1 less than the number of vertices
   if(SltEdgNum != SltVtxNum-1)   
   {   
     cout << "\n The current solution is unfeasible: the number of edges does not match the number of the vertices!"<<endl;
     getchar();
     return false;                     
   }     
   //Check if the budget constraint is satisfied     
   if(GetSolutionBudget()> Max_Budget)
   {
     cout<<"\n The current solution is unfeasible: the budget constraint is not satisfied!" <<endl;
     getchar();
     return false;                           
   }
   
   for(int i=1;i<Max_Vtx;i++)  
     if(IfVtxOnSolution[i])  
     {
        int CurVtxDepth=GetVtxDepth(i);
        if(CurVtxDepth>=MaxHop) 
        {
           cout<<"\n The current solution is unfeasible: the Hop-Constraint is not satisfied!" <<endl;
           getchar();
           return false;                        
        } 
        if(CurVtxDepth==Null)
        {
           cout << "\n The current solution is unfeasible: the graph is not connected! unconnected node:"<<i+1<<endl;
           getchar();
           return false;
        }                                      
     }  
 
   if(IfDisplayDetail)
   {  
      cout<<" The total revenue of the current solution is:"<< GetSolutionRevenue()<<endl;
      cout<<" The total cost of the current solution is:"<< GetSolutionBudget()/100<<endl;  
      cout<<" The current solution is feasible!"<<endl;    
   }
   //cout<<" The current solution is feasible!"<<endl;
   return true;    
}//End CheckSolutionFeasible

//Calcult the cost of the shortest path from the root to every node, with hop limited to be less than MaxHop
void GetShortestPathWithHop(int LimitedHop)
{
   if(LimitedHop > MaxHop)
   {
      cout <<"\n GetShortestPathWithHop(): Error! the limited hop is larger than the max hop"<<endl;      
      return;
   }
     
   ShortestPath[0][0]=0;
   for(int i=1;i<Max_Vtx;i++)
     ShortestPath[i][0]=InfiniteEdgeCost;
         
   for(int i=1;i<LimitedHop;i++) 
      for(int j=0;j<Max_Vtx;j++)   
      {  
         ShortestPath[j][i]=ShortestPath[j][i-1];     
         for(int k=0;k<Max_Vtx;k++)
         {
            if(ShortestPath[k][i-1]+EdgeBudget[k][j]<ShortestPath[j][i])
              ShortestPath[j][i]=ShortestPath[k][i-1]+EdgeBudget[k][j];                                             
         }
      } 
}//End GetShortestPathCostWithHop()

//Calcult the cost of the shortest path from the root to every node, with hop limited to be less than MaxHop
void GetCompleteShortestPathWithHop(int LimitedHop)
{     
   starting_time = (double)clock();   
   if(IfDisplayDetail)
     cout<<"Begin to Get Complete Shortest path"<<endl;    
   if(LimitedHop > MaxHop)
   {
      cout <<"\nGetCompleteShortestPath(): Error! the limited hop is larger than the max hop"<<endl;
      getchar();
      return;
   }
        
   CompleteShortestPath[0][0][0]=0;
   for(int i=0;i<Max_Vtx;i++)
     for(int j=0;j<Max_Vtx;j++)
       if(i==j)
         CompleteShortestPath[i][j][0]=0;
       else
         CompleteShortestPath[i][j][0]=InfiniteEdgeCost;   
  
   for(int i=0; i<Max_Vtx; i++)
   {          
     //if(VtxRevenue[i]>0)
     {
       SwapVertex(0,i);
   
       GetShortestPathWithHop(LimitedHop);    
       for(int j=0;j<Max_Vtx;j++)
         for(int k=0; k<MaxHop; k++)
         {
           if(j==0)
             CompleteShortestPath[i][j][k]=ShortestPath[i][k];
           else if(j==i)
             CompleteShortestPath[i][j][k]=ShortestPath[0][k];
           else
             CompleteShortestPath[i][j][k]=ShortestPath[j][k];          
         }
           
         SwapVertex(0,i); 
       }
    }//End for     
 
   finishing_time = (double)clock(); 
   PrepareTime=finishing_time-starting_time; 
   cout<<"Finish to Get Complete Shortest path. Used time:"<<PrepareTime<<endl;  
}//End GetShortestPathCostWithHop()

void SaveShortestPath()
{
   for(int i=0; i<Max_Vtx; i++)
     for(int j=0; j<MaxPathLength; j++)
     TempShortestPath[i][j]=ShortestPath[i][j];   
}//End SaveShortestPath()

void RestoreShortestPath()
{
   for(int i=0; i<Max_Vtx; i++)
     for(int j=0; j<MaxPathLength; j++)
     ShortestPath[i][j]=TempShortestPath[i][j];
}//End RestoreShortestPath()

void SaveShortestPathToCurSlt()
{
   for(int i=0; i<Max_Vtx; i++)    
     TempShortestPathToCurSlt[i]=ShortestPathToCurSlt[i];   
}//End SaveShortestPathToCurSlt()

void RestoreShortestPathToCurSlt()
{
   for(int i=0; i<Max_Vtx; i++)     
     ShortestPathToCurSlt[i]=TempShortestPathToCurSlt[i];
}//End RestoreShortestPathToCurSlt()

void GetShortestPathToCurSlt()
{
   UpdateSolution();
   for(int i=0; i<Max_Vtx; i++)
   {
      ShortestPathToCurSlt[i]=InfiniteEdgeCost;
      ShortestPathConnectedVtx[i]=0;
      if(IfVtxOnSolution[i])
      {
         ShortestPathToCurSlt[i]=0;         
         continue;
      }
      
      if(VtxRevenue[i]==0)
      {
         ShortestPathToCurSlt[i]=Null;   
         continue;              
      }      
      
      for(int j=1; j<=VtxOnSolution[0]; j++)
      {        
        if(ShortestPathToCurSlt[i] > ShortestPathMatrix[i][VtxOnSolution[j]]) 
        { 
          ShortestPathToCurSlt[i]=ShortestPathMatrix[i][VtxOnSolution[j]]; 
          ShortestPathConnectedVtx[i]=VtxOnSolution[j];
        }     
      }//End for            
   }//End for      
}//End GetShortestPathToCurSlt

void SortFitnissVtx()
{
   VtxWithFitNum=0;
   for(int i=0; i<Max_Vtx; i++)
   {
      if( VtxRevenue[i]==0 || VtxStatus[i]==2 || IfVtxOnSolution[i] || ShortestPath[i][MaxHop-1] >= InfiniteEdgeCost)
      {        
        VtxFitniss[i]=0;
        continue;
      }
      else
      {
        CandidateVtx[VtxWithFitNum++]=i;        
        int PathCost =ShortestPathToCurSlt[i];   
        VtxFitniss[i] = (float)(VtxRevenue[i]*VtxRevenue[i]*VtxRevenue[i])/PathCost; 
       }   
   } 
  
  int TempVtx;
  for(int i=0; i<VtxWithFitNum-1;i++)
    for(int j=i+1; j<VtxWithFitNum;j++)
    {
       if(VtxFitniss[CandidateVtx[i]]<VtxFitniss[CandidateVtx[j]])
       {
          TempVtx=CandidateVtx[i];
          CandidateVtx[i]=CandidateVtx[j];
          CandidateVtx[j]=TempVtx;                                          
       }        
    } 
}//End SortFitnissVtx()

//Init an empty solution which only containes the root and does not contain any edge
void InitEmptySolution()
{
   for(int i=0; i<Max_Vtx; i++)
   {
      FinalSolution[i]=Null;     
   } 
   
   VtxStatus[0] = 1;  
   for(int i=1; i<Max_Vtx; i++)
     VtxStatus[i] = 0;
  
  GetShortestPathToCurSlt();       
}//End InitEmptySolution

void UpdShortestPathToCurSltAftDel(int DeleteVtxNum)
{     
   ShortestPathToCurSlt[DeleteVtxNum]=Null;
   for(int i=0; i<Max_Vtx; i++)
   {     
      if(i==DeleteVtxNum || IfVtxOnSolution[i])     
        continue;      
      
      if(ShortestPathToCurSlt[i] == ShortestPathMatrix[i][DeleteVtxNum])
      {
         ShortestPathToCurSlt[i]=InfiniteEdgeCost;
         for(int j=0; j<Max_Vtx; j++)
         {         
           if(j!=DeleteVtxNum && IfVtxOnSolution[j] && ShortestPathToCurSlt[i] > ShortestPathMatrix[i][j])  
           {
              
              ShortestPathToCurSlt[i]=ShortestPathMatrix[i][j];              
           }  
         }//End for 
      }      
   }     
}//End UpdShortestPathToCurSltAftDel()

//Find the last node iteratively of the hortest path from the root to the current node (with less than LimitedHop hops) 
bool FindLastNode(int EndNode, int BeginNode)
{    
   RecordPath[PathLength++] = EndNode;     
   if(EndNode==BeginNode)    
     return true;      
   
   int LastNode;  
   for(int i=1;i<=AdjacentVtx[EndNode][0];i++)     
   { 
      LastNode = AdjacentVtx[EndNode][i];                 
      if(LastNode!=EndNode && ShortestPathMatrix[BeginNode][LastNode]+TempEdgeBudget[LastNode][EndNode]==ShortestPathMatrix[BeginNode][EndNode]) 
      {                                       
         FindLastNode(LastNode, BeginNode);  
           break;      
      }  
   }  
     
   if(PathLength>1 && RecordPath[PathLength-1]==BeginNode)
     return true;
   else     
     return false;       
}//End FindLastNode

//Get the shortest path from the root to a node (with less than LimitedHop hops)
bool GetRealPathBetweenTwoVtx(int EndNode, int BeginNode)
{      
   PathLength =0;
   for(int i=0; i<MaxPathLength; i++) 
     RecordPath[i]=-1;
    
   if(FindLastNode(EndNode, BeginNode)) 
     return true;             
   else 
   {
     cout<<"\n Fail to get real path from node "<<SimplyVtxIndex[BeginNode]+1<< " to node "<<SimplyVtxIndex[EndNode]+1<<endl; 
     return false;  
   }
}//End GetRealPathBetweenTwoVtx

//Evalute the fitness function of a path
float EvalutePath(int EndNode)
{
   if( VtxRevenue[EndNode]==0 ||VtxStatus[EndNode]==2 || IfVtxOnSolution[EndNode] || ShortestPath[EndNode][MaxHop-1] >= InfiniteEdgeCost)
     return 0;  

   int PathCost =ShortestPathToCurSlt[EndNode];  
   if(PathCost >0 && PathCost + CurSltBudget > Max_Budget)   
   {  
      return 0;
   }
   
   PathCost =ShortestPathMatrix[0][EndNode];
   int BeginNode=Null;
   //if(!IfVtxOnSolution[BeginNode])
   {       
      for(int i=0; i<Max_Vtx; i++)
      {
        if(i!=EndNode && IfVtxOnSolution[i] && ShortestPathMatrix[i][EndNode] <= PathCost)  
        {
           BeginNode=i; 
           PathCost=ShortestPathMatrix[i][EndNode];
        }     
      } 
   }  
 
   if(!GetRealPathBetweenTwoVtx(EndNode, BeginNode))
     return 0;
  
   int EdgeNumBtwTwoVtx = PathLength-1; 
   if(GetVtxDepth(BeginNode)+ EdgeNumBtwTwoVtx > MaxHop-1)
     return 0;        
  
   if(PathCost >0 && PathCost + CurSltBudget <= Max_Budget)
   //if(PathCost >0 && PathCost + GetSolutionBudget() <= Max_Budget)
   {  
      float CurNodeScore = (float)(VtxRevenue[EndNode]*VtxRevenue[EndNode]*VtxRevenue[EndNode])/PathCost; 
      return CurNodeScore;
   }
   else if(PathCost <=0)
   {
      cout <<"\n Error! The cost from the root to vertex "<<SimplyVtxIndex[EndNode]+1 << " equals to 0, please check it!" << endl;    
      getchar();
      return 0;
   }
   else    
      return 0; 
}//End EvalutePath()

bool AddOnePathToSolution(int EndNode)
{   
   if(IfVtxOnSolution[EndNode])
   {
     cout<<"Error! AddOnePathToSolution: node "<<SimplyVtxIndex[EndNode]+1<<" is already on the solution"<<endl;
     getchar();
     return false;                            
   }
   
   int ShortestPathCost = ShortestPathMatrix[0][EndNode];  
   int BeginNode=0;
   int EdgeNumBtwTwoVtx;
 
   for(int i=1; i<Max_Vtx; i++)      
   {
      if(IfVtxOnSolution[i] && ShortestPathMatrix[i][EndNode] <= ShortestPathCost)
      {           
         ShortestPathCost= ShortestPathMatrix[i][EndNode];
         BeginNode = i;                          
      }        
   }
  
   if(!GetRealPathBetweenTwoVtx(EndNode, BeginNode))
     return false; 
  
   EdgeNumBtwTwoVtx = PathLength-1; 
   if(GetVtxDepth(BeginNode)+ EdgeNumBtwTwoVtx > MaxHop-1)
   {
      cout<<"The distance between two vertice is too large:"<<SimplyVtxIndex[BeginNode]+1<<" "<<SimplyVtxIndex[EndNode]+1<<" "<<GetVtxDepth(BeginNode)<<" "<<EdgeNumBtwTwoVtx<<endl;
      return false;  
   }
  
   int TempSltBudget = CurSltBudget;  
   for(int i=0; i<PathLength-1; i++)
   {
      int CurrentNode =RecordPath[i];
      int LastNode=RecordPath[i+1];
      if(FinalSolution[CurrentNode]!=LastNode)   
        TempSltBudget += TempEdgeBudget[LastNode][CurrentNode];  
   }
   if(TempSltBudget > Max_Budget)
   {
      cout<<"The budget is too large:"<<ShortestPathCost<<" "<<TempSltBudget-CurSltBudget<<" "<<ShortestPathToCurSlt[EndNode]<<endl;
      getchar();
      return false;
   }
        
    for(int i=0; i<PathLength-1; i++)
    {
       int CurrentNode =RecordPath[i];
       int LastNode=RecordPath[i+1];
       if(!IfVtxOnSolution[CurrentNode])
       {
          FinalSolution[CurrentNode]=LastNode;          
          IfVtxOnSolution[CurrentNode]=true;          
          IfLeafVtx[CurrentNode]=false;
          IfLeafVtx[LastNode]=false; 
          CurSltBudget += EdgeBudget[LastNode][CurrentNode];
          SltVtxNum++;
          SltEdgNum++;   
      }    
    } 
    IfLeafVtx[RecordPath[0]]=true;  
 
    if(GetVtxDepth(EndNode)>=MaxHop)
    {
      cout <<"The depth of node "<<SimplyVtxIndex[EndNode]+1<<" is "<<GetVtxDepth(EndNode)<<endl;   
      getchar();     
    }
    return true; 
}//End AddOnePathToSolution()

int GetBestVtxWithoutHop()
{
   int BestVtxNum=Null;
   int LastValidVtx=Null; 
   for(int i=0; i<VtxWithFitNum; i++)
   {
       if(EvalutePath(CandidateVtx[i])>0)
       {
          LastValidVtx=CandidateVtx[i];
          if(RandomInt(100)<=PrabToAddBestVtx)
          {
             BestVtxNum=CandidateVtx[i];   
             break;
          }                               
       }    
   }
   
  if(BestVtxNum==Null)
     BestVtxNum=LastValidVtx;         
   
   return BestVtxNum;     
}//End GetBestVtxWithoutHop()

bool AddBestPathWithoutHop()
{
   int BestNode=GetBestVtxWithoutHop();
   if(BestNode!=Null && AddOnePathToSolution(BestNode))  
     return true;
     
   return false;     
}//End AddBestPathWithoutHop()

//Find the last node iteratively of the hortest path from the root to the current node (with less than LimitedHop hops) 
bool IterFindLastNode(int EndNode, int BeginNode, int LimitedHop)
{
   if(LimitedHop > MaxHop)
   {
      cout <<"\n FindLastNode(): Error! the limited hop is larger than the max hop"<<endl;
      return false;
   }
    
   RecordPath[PathLength++] = EndNode;     
   if(EndNode==BeginNode)    
     return true;      
   
   int LastNode;  
   for(int i=1;i<=AdjacentVtx[EndNode][0];i++)
   { 
      LastNode = AdjacentVtx[EndNode][i];                
      if(LastNode!=EndNode && CompleteShortestPath[BeginNode][LastNode][LimitedHop-2]+EdgeBudget[LastNode][EndNode] == CompleteShortestPath[BeginNode][EndNode][LimitedHop-1]) 
      {                                       
          IterFindLastNode(LastNode, BeginNode, LimitedHop-1);  
          break;      
       }  
   }
     
   if(PathLength>1 && RecordPath[PathLength-1]==BeginNode)
      return true;
   else     
      return false;       
}//End FindLastNode

//Get the shortest path from the root to a node (with less than LimitedHop hops)
bool GetRealPathBetweenTwoVtx(int EndNode, int BeginNode, int LimitedHop)
{ 
   if(LimitedHop > MaxHop)
   {
     cout <<"\n GetRealPath(): Error! the limited hop is larger than the max hop"<<endl;
     getchar();
     return false;
   }
     
   PathLength =0;
   for(int i=0; i<MaxPathLength; i++) 
     RecordPath[i]=-1;
    
   if(IterFindLastNode(EndNode, BeginNode, LimitedHop)) 
     return true;             
   else 
   {
     cout<<"\n Fail to get real path from the root to node "<<EndNode+1<<endl; 
     return false;  
   }
}//End GetRealPath

//Evalute the fitness function of a path
float EvalutePath(int EndNode, int LimitedHop)
{
   if(VtxRevenue[EndNode]==0 ||VtxStatus[EndNode]==2 || IfVtxOnSolution[EndNode] || ShortestPath[EndNode][LimitedHop-1] >= InfiniteEdgeCost)
     return 0;  

   int PathCost=CompleteShortestPath[EndNode][0][LimitedHop-1];
   int BeginNode=Null;
   for(int i=0; i<Max_Vtx; i++)
   {
      if(i!=EndNode && IfVtxOnSolution[i] && CompleteShortestPath[EndNode][i][LimitedHop-1-GetVtxDepth(i)] <= PathCost) 
      {
         BeginNode=i; 
         PathCost=CompleteShortestPath[EndNode][i][LimitedHop-1-GetVtxDepth(i)];             
      }        
   } 
  
   if(PathCost >0 && PathCost + CurSltBudget <= Max_Budget)   
   {  
      float CurNodeScore = (float)(VtxRevenue[EndNode]*VtxRevenue[EndNode]*VtxRevenue[EndNode])/PathCost; 
      return CurNodeScore;
   }
   else if(PathCost<=0)
   {
      cout <<"\n Error! The cost from the root to vertex "<<SimplyVtxIndex[EndNode]+1 << " equals to 0, please check it!" << endl;    
      getchar();
      return 0;
   }
   else  
      return 0;         
}//End EvalutePath()

int GetBestVtxWithHop(int LimitedHop)
{
   int BestVtxNum=Null;
   int LastValidVtx=Null; 
   for(int i=0; i<VtxWithFitNum; i++)
   {
       if(EvalutePath(CandidateVtx[i],LimitedHop)>0)
       {
          LastValidVtx=CandidateVtx[i];
          if(RandomInt(100)<=PrabToAddBestVtx)
          {
             BestVtxNum=CandidateVtx[i];   
             break;
          }                               
       }    
   }//End For
   
  if(BestVtxNum==Null)
     BestVtxNum=LastValidVtx;         
   
   return BestVtxNum;        
}//End GetBestVtxWithHop

int AddOnePathToSolution(int BestNode, int LimitedHop)
{
   if(IfVtxOnSolution[BestNode])
   {
     cout<<"Error! AddOnePathToSolution: node "<<SimplyVtxIndex[BestNode]+1<<" is already on the solution"<<endl;
     getchar();
     return false;                            
   }
   
   int ShortestPathCost = CompleteShortestPath[BestNode][0][LimitedHop-1];  
   int BeginNode=0;
   int EdgeNumBtwTwoVtx;
 
   for(int i=1; i<Max_Vtx; i++)      
   {
      if(i!=BestNode && IfVtxOnSolution[i] && CompleteShortestPath[BestNode][i][LimitedHop-1-GetVtxDepth(i)] <= ShortestPathCost)
      {           
         ShortestPathCost= CompleteShortestPath[BestNode][i][LimitedHop-1-GetVtxDepth(i)];
         BeginNode = i;                          
      }        
   }

   if(!GetRealPathBetweenTwoVtx(BestNode, BeginNode, LimitedHop-GetVtxDepth(BeginNode)))
     return false; 
  
   EdgeNumBtwTwoVtx = PathLength-1;   
   int TempSltBudget = CurSltBudget;
   
   bool IfLoopExist=false;  
   for(int i=0; i<PathLength-1; i++)
   {
      int CurrentNode =RecordPath[i];
      if(IfVtxOnSolution[CurrentNode])
      {
         if(IfDisplayDetail)
           cout<<"loop exist"<<endl;
         IfLoopExist=true;
      }
      int LastNode=RecordPath[i+1];
      if(FinalSolution[CurrentNode]!=LastNode)   
        TempSltBudget += TempEdgeBudget[LastNode][CurrentNode];  
   }
   if(TempSltBudget>Max_Budget)
   {
      cout<<"The budget is too large:"<<ShortestPathCost<<" "<<TempSltBudget-CurSltBudget<<" "<<ShortestPathToCurSlt[BestNode]<<endl;
      getchar();
      return false;
   }
 
    if(!IfLoopExist)
    {
       for(int i=0; i<PathLength-1; i++)
       {
         int CurrentNode =RecordPath[i];
         int LastNode=RecordPath[i+1];
         if(!IfVtxOnSolution[CurrentNode])
         {
           FinalSolution[CurrentNode]=LastNode;          
           IfVtxOnSolution[CurrentNode]=true;          
           IfLeafVtx[CurrentNode]=false;
           IfLeafVtx[LastNode]=false; 
           CurSltBudget += EdgeBudget[LastNode][CurrentNode];
           SltVtxNum++;
           SltEdgNum++;   
         }    
       } 
       IfLeafVtx[RecordPath[0]]=true;                     
    }
    else
    {
       for(int i=0; i<PathLength-1; i++)
       {
         int CurrentNode =RecordPath[i];
         int LastNode=RecordPath[i+1];      
         FinalSolution[CurrentNode]=LastNode; 
       }  
       
       UpdateSolution();         
    } 
 
    if(GetVtxDepth(BestNode)>=MaxHop)
    {
      cout <<"The depth of node "<<SimplyVtxIndex[BestNode]+1<<" is "<<GetVtxDepth(BestNode)<<endl;   
      getchar();     
    }
    return true;                                
}//End AddOnePathToSolution

bool AddBestPathWithHop(int LimitedHop)
{     
   int BestNode=GetBestVtxWithHop(LimitedHop);
   if(BestNode!=Null && AddOnePathToSolution(BestNode, LimitedHop))  
     return true;
     
   return false;  
}//End AddBestPathWithHop

//Add the best path to the currrent solution
bool AddBestPathToSolution()
{  
   if(AddBestPathWithoutHop()) 
     return true;  

   if(CurSltBudget>=Max_Budget-100)
     return false;
   
   int UncollectedVtxNum=0;
   for(int i=0; i<Max_Vtx; i++)
   if( VtxRevenue[i]>0 && !IfVtxOnSolution[i] && VtxStatus[i]!=2 && ShortestPath[i][MaxHop-1] < InfiniteEdgeCost)
     UncollectedVtxNum++;  
    
   if(UncollectedVtxNum==0)
      return false;    

   if(AddBestPathWithHop(MaxHop))
     return true; 
      
   return false;
}

//Construct an initial solution by iteratively add the best path into the current solution
bool ConstructInitialTree()
{       
   for(int i=0; i<Max_Vtx; i++)
     for(int j=0; j<Max_Vtx; j++)
     EdgeBudget[i][j] = TempEdgeBudget[i][j];     
     
   RestoreShortestPath();  
   
   InitEmptySolution(); 
   GetShortestPathToCurSlt();     
   SortFitnissVtx();  
   while(AddBestPathToSolution())
   {
     GetShortestPathToCurSlt();     
     SortFitnissVtx();                              
   };    
   UpdateSolution();
   if(CheckSolutionFeasible())   
     return true;
 
   cout<<"\n Fail to construct an initial solution"<<endl;
   getchar();
   return false;  
}//End ConstructInitialTree

float GetFitnessOfCurSolution()
{      
   int CurSltBudget = GetSolutionBudget();    
   float CurSltRevenue=(float)GetSolutionRevenue();
   if(Alpha==Null)
   {
      if(CurSltBudget <= Max_Budget)
        return CurSltRevenue*Max_Budget-CurSltBudget;      
      else
        return Null;  
   } 
   else if(Alpha>=0)
   {
      if(CurSltBudget<=0 || CurSltBudget > Max_Budget)
        return Null;
      float CurSltfit=1/(float)CurSltBudget;
      for(int i=0; i<Alpha; i++)
        CurSltfit*=(float)CurSltRevenue/100;           
      return CurSltfit;      
   }
   else
     return Null; 
}//GetFitnessOfCurSolution()

void SaveCurSltToTempArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     TempSolutionArray[i]=FinalSolution[i];
}//End SaveCurSltToTempArray()

void RestoreSolutionFromTempArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     FinalSolution[i]=TempSolutionArray[i];  
       
   UpdateSolution();
}//RestoreSolutionFromTempArray()

void SaveCurSltToBestArray()
{ 
   for(int i=0; i<Max_Vtx; i++)     
     BestSolutionArray[i]=FinalSolution[i];     
}//End SaveCurSltToBestArray()

void RestoreSolutionFromBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     FinalSolution[i]=BestSolutionArray[i];   
 
   UpdateSolution(); 
}//End RestoreSolutionFromBestArray()

void SaveCurSltToCurIterBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     CurIterBestSolutionArray[i]=FinalSolution[i];   
}//End SaveCurSltToCurIterBestArray()

void RestoreSolutionFromCurIterBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     FinalSolution[i]=CurIterBestSolutionArray[i];   
 
   UpdateSolution(); 
}//End RestoreSolutionFromCurIterBestArray()

void SaveCurSltToBLSBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     BLSBestSolutionArray[i]=FinalSolution[i];  
}//End SaveCurSltToGlobalBestArray()

void RestoreSolutionFromBLSBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     FinalSolution[i]=BLSBestSolutionArray[i];   
 
   UpdateSolution(); 
}//End RestoreSolutionFromGlobalBestArray()

void SaveCurSltToGlobalBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     GlobalBestSolutionArray[i]=FinalSolution[i];  
}//End SaveCurSltToGlobalBestArray()

void RestoreSolutionFromGlobalBestArray()
{
   for(int i=0; i<Max_Vtx; i++)     
     FinalSolution[i]=GlobalBestSolutionArray[i];   
 
   UpdateSolution(); 
}//End RestoreSolutionFromGlobalBestArray()

void StoreCurrentSolution(int StoredIndex)
{
   if(StoredIndex<0 || StoredIndex >= MaxStoredSolutionNum)    
     return;
   
   for(int i=0; i<Max_Vtx;i++)   
     StoredSolution[StoredIndex][i]=FinalSolution[i];     
}//End StoreCurrentSolution

void RestoreSolution(int StoredIndex)
{
   if(StoredIndex<0 || StoredIndex >= MaxStoredSolutionNum) 
     return;                
     
   for(int i=0; i<Max_Vtx;i++)
   {
     FinalSolution[i]=StoredSolution[StoredIndex][i];      
   } 
   UpdateSolution();
}//End RestoreSolution

int CalDistBtwCurSltAndStoredSlt(int StoredIndex)
{
   if(StoredIndex<0 || StoredIndex >= MaxStoredSolutionNum)    
      return Null;   
   
   int Distance=0;
   for(int i=0; i<Max_Vtx; i++)  
      if(FinalSolution[i]!=StoredSolution[StoredIndex][i])
        Distance++;
        
   return Distance;    
}//End CalDistBtwTwoSolution

int CalDistBtwTwoSolution(int FirstStoredIndex, int SecondStoredIndex)
{
   if(FirstStoredIndex<0 || FirstStoredIndex >= MaxStoredSolutionNum || SecondStoredIndex<0 || SecondStoredIndex >= MaxStoredSolutionNum)    
     return Null;  
   
   int Distance=0;
   for(int i=0; i<Max_Vtx; i++)  
      if(StoredSolution[FirstStoredIndex][i]!=StoredSolution[SecondStoredIndex][i])
        Distance++;
        
   return Distance;    
}//End CalDistBtwTwoSolution

int GetVtxOutputNum(int CurVtx)
{
   int OutputNum =0;  
   
   for(int i=1; i<=VtxOnSolution[0]; i++)    
     if(CurVtx!=VtxOnSolution[i] && FinalSolution[VtxOnSolution[i]]==CurVtx)   
       OutputNum++;    

   return OutputNum;      
}

bool DeleteLeafVtx(int LeafVtxNum)
{
   if(LeafVtxNum==0 || !IfLeafVtx[LeafVtxNum]) 
     return false;                      

   int CurrentNode = LeafVtxNum;
   int ParentNode = FinalSolution[CurrentNode];

   do
   {
      FinalSolution[CurrentNode]=Null; 
      IfVtxOnSolution[CurrentNode]=false;       
      IfVtxOnSolution[CurrentNode]=false;
      IfLeafVtx[CurrentNode]=false;
      if(GetVtxOutputNum(ParentNode)>0)
        IfLeafVtx[ParentNode]=false;
      else
        IfLeafVtx[ParentNode]=true;
        
      SltEdgNum--;
      SltVtxNum--;
      CurSltBudget -= EdgeBudget[ParentNode][CurrentNode];
              
      CurrentNode = ParentNode;
      ParentNode = FinalSolution[CurrentNode];     
   }while(CurrentNode!=-1 && VtxRevenue[CurrentNode]==0 && GetVtxOutputNum(CurrentNode)==0);

   return true;
}//End DeleteLeafVtx()

bool AutoDeleteOneLeaf()
{   
   int BeginVtxNum=RandomInt(Max_Vtx);
   int DeleteLeafNum=BeginVtxNum;   
   for(int i=0;i<Max_Vtx;i++)
   {
      DeleteLeafNum=(BeginVtxNum+i)%Max_Vtx;      
      if(IfLeafVtx[DeleteLeafNum] && RandomInt(RecordSltNum)>VtxOccurTimes[DeleteLeafNum]) 
        break;    
   }
   
   if(IfLeafVtx[DeleteLeafNum]) 
   { 
      DeleteLeafVtx(DeleteLeafNum);   
      return true;
   }
   
   int WorstLeafVtx=0;
   for(int i=1;i<Max_Vtx;i++)
   {      
      if(IfLeafVtx[i] && VtxOccurTimes[i]<VtxOccurTimes[WorstLeafVtx])
        WorstLeafVtx=i;     
   }   
   
   if(WorstLeafVtx!=0 && IfLeafVtx[WorstLeafVtx])
   {
      DeleteLeafVtx(WorstLeafVtx);   
      return true;
   }
   else
      return false;   
}

void ReconstructFromBackBone()
{
  while(AddBestPathToSolution())
    ;  
}//End ReconstructFromBackBone()

void AddAfterDeleteLeaf(int DeleteLeafVtxNum)
{ 
   DeleteLeafVtx(DeleteLeafVtxNum);
   ReconstructFromBackBone(); 
}//End AddafterDeleteLeaf

void AddAfterDeleteTwoLeaf(int FirstLeafVtxNum, int SecondLeafVtxNum)
{
   DeleteLeafVtx(FirstLeafVtxNum);  
   DeleteLeafVtx(SecondLeafVtxNum); 
   ReconstructFromBackBone(); 
}//End AddafterDeleteLeaf

bool AddAfterDeleteSeveralVtx(int DeleteVtxNum)
{ 
   if(DeleteVtxNum==0)
     return true;

   int cnt=0;
   bool IfDeleteVtx=false;
   while(cnt++<DeleteVtxNum && AutoDeleteOneLeaf())  
     IfDeleteVtx=true;

   if(IfDeleteVtx)
     ReconstructFromBackBone();  
     
   return IfDeleteVtx;      
}//End AddAfterDeleteSeveralVtx()

void GenerateNS(int FisrtDeleteVtx, int DeleteVtxNum)
{ 
   DeleteLeafVtx(FisrtDeleteVtx);  
   AddAfterDeleteSeveralVtx(DeleteVtxNum);  
}//End GenerateNS

void NeighborhoodSearch(int MaxVtxNumToDelete)
{   
   float CurSltFitness; 
   float BestFitness; 
   bool HasFoundBetterSolution;
   int FirstIndex, FirstLeaf, SecondIndex, SecondLeaf;
   bool  IfTestDelete[MaxVtxNum];
   for(int i=0; i<MaxVtxNum;i++)
     IfTestDelete[i]=false;     
     
   SaveCurSltToBestArray(); 
   SaveCurSltToCurIterBestArray();    
     
   BestFitness = GetFitnessOfCurSolution();  
   do
   {               
      BeginNewIter: HasFoundBetterSolution=false;   
      RestoreSolutionFromCurIterBestArray(); 
      SaveCurSltToBestArray();         
      RestoreSolutionFromBestArray();        
      GetShortestPathToCurSlt();          
      SaveShortestPathToCurSlt();     
      SortFitnissVtx();    
      for(int i=0; i<Max_Vtx;i++)
      {
         if(IfLeafVtx[i])
           IfTestDelete[i]=true;
         else
           IfTestDelete[i]=false;
      }  
      
      FirstIndex=RandomInt(Max_Vtx);                   
      for(int i=0; i<Max_Vtx; i++)
      {
         FirstLeaf = (FirstIndex+i)%Max_Vtx;    
         if(!IfTestDelete[FirstLeaf])  
           continue;      
         
         for(int m=0; m<MaxVtxNumToDelete;m++) 
         {
            RestoreSolutionFromBestArray(); 
            RestoreShortestPathToCurSlt();  

            IfDisplayDetail=false;             
            GenerateNS(FirstLeaf,m);    
 
            CurSltFitness = GetFitnessOfCurSolution();          
            if(CurSltFitness > BestFitness)    
            {              
              SaveCurSltToCurIterBestArray();          
              BestFitness = CurSltFitness; 
              HasFoundBetterSolution=true;
            }//End if() 
         }//End for(int m=0; m<MaxVtxToDelete,m++)    
      }//End for(int i=0; i<Max_Vtx; i++) 
   }while(HasFoundBetterSolution);  
     
   RestoreSolutionFromBestArray();  
}//End NeighborhoodSearch()

void Perturb(int PerturbRank)
{
  UpdateSolution();
  GetShortestPathToCurSlt();          
  SaveShortestPathToCurSlt();     
  SortFitnissVtx();
  int DeleteVtxNum=PerturbRank;
  AddAfterDeleteSeveralVtx(DeleteVtxNum);
}//End Perturb()

int GetCollectedProfitVtxNum()
{
   int CollectedProfitVtxNum=0;
   for(int i=1; i<Max_Vtx; i++)
   {
     if(VtxRevenue[i]>0 && IfVtxOnSolution[i])  
       CollectedProfitVtxNum++;      
   }
   return CollectedProfitVtxNum;
}// End GetCollectedProfitVtxNum()

void RecordVtxOccurTimes()
{
   for(int i=0; i<Max_Vtx; i++)
     VtxOccurTimes[i]=0;     
   
   int AvgRevenue=0;
   int LimitedRevenue=0;
   int TempRecordSltNum=0;
   for(int i=0; i<10; i++)
   { 
      ConstructInitialTree();   
      NeighborhoodSearch(1);       
      AvgRevenue+=GetSolutionRevenue();  
   }   
   LimitedRevenue=AvgRevenue/10; 
 
   RecordSltNum=0;
   int TempTotalRevenue=0;
   AvgDistance=0;   
   int TempCnt=0;
   while(RecordSltNum<MaxRecordSolutionNum)
   {     
      IfDisplayDetail=false;     
      ConstructInitialTree();  
      NeighborhoodSearch(1);
      UpdateSolution();    
      
      TempCnt++;
      if(TempCnt%2==0)
        StoreCurrentSolution(MaxStoredSolutionNum-1); 
      else
        AvgDistance+=CalDistBtwCurSltAndStoredSlt(MaxStoredSolutionNum-1);
    
      if(GetSolutionRevenue()>=LimitedRevenue)
      { 
         RecordSltNum++;  
         TempTotalRevenue+=GetSolutionRevenue();                      
         for(int k=0; k<Max_Vtx; k++)     
           if(IfVtxOnSolution[k]) 
             VtxOccurTimes[k]++;
      }
   }  
   
   AvgDistance=2*AvgDistance/RecordSltNum;      
   return;
}

int GetAverageDistance(int MaxSltNum)
{
   if(MaxSltNum==0)
      return Null;

   int  SumDis=0;
   int TempTotalRevenue=0;
   int TempCnt=0;
   for(int i=0; i<MaxSltNum; i++)
   {  
      ConstructInitialTree();   
      NeighborhoodSearch(MaxNS);    
      StoreCurrentSolution(MaxStoredSolutionNum-1);    
      ConstructInitialTree();   
      NeighborhoodSearch(MaxNS); 
      SumDis+=CalDistBtwCurSltAndStoredSlt(MaxStoredSolutionNum-1);  
   }
   return SumDis/MaxSltNum;
}

//Breakout local search procedure
void BLS()
{    
   double BLSBeginTime = (double)clock();   
   int Distance=0;   
   int MinPerturbRank=1;   
   int MaxPerturbRank;   
   int PerturbRank=MinPerturbRank;  
   int BestRevenue;
                
   ConstructInitialTree();      
   NeighborhoodSearch(MaxNS);    
   SaveCurSltToBLSBestArray();
   
   BestRevenue=GetSolutionRevenue();   
   MaxPerturbRank=GetCollectedProfitVtxNum();
    
   float BestSltFit=GetFitnessOfCurSolution();
   float NewSltFit;
   int NoImprovePerturbation=0;  

   while(GetSolutionRevenue()<ExpectedRevenue && NoImprovePerturbation<MaxNoImpPerturb && (double)(clock())-BLSBeginTime < MaxUsedTime)
   {       
      StoreCurrentSolution(MaxStoredSolutionNum-1);        
   
      if(!AddAfterDeleteSeveralVtx(PerturbRank))
        break;
      NeighborhoodSearch(MaxNS);
                  
      Distance=CalDistBtwCurSltAndStoredSlt(MaxStoredSolutionNum-1);
      MaxPerturbRank=GetCollectedProfitVtxNum();   
 
      if(Distance<=LimitedDis*AvgDistance)
      {   
         PerturbRank+=1;           
         if(PerturbRank>MaxPerturbRank) 
           PerturbRank=MaxPerturbRank;    
      }
      else
      { 
         PerturbRank-=1;
         if(PerturbRank<MinPerturbRank)       
           PerturbRank=MinPerturbRank; 
      }//End if()        
      
      NewSltFit=GetFitnessOfCurSolution();          
      if(NewSltFit>BestSltFit)
      {
         SaveCurSltToBLSBestArray();
         BestSltFit = NewSltFit;
         BestRevenue=GetSolutionRevenue();
         NoImprovePerturbation=0;         
      } 
      else
        NoImprovePerturbation++;             
   }//End while(NoImprovePerturbation <20...)  

   RestoreSolutionFromBLSBestArray();   
}//End BLS()

void Multi_BLS()
{   
   int MultiStartTimes=0; 
   int GlobalBestRevenue=0;  
   int CurSltRevenue;
   int TotalRevenue=0;     
     
   while(MultiStartTimes++< RunBLSTimes)	
   {      
      starting_time = (double)clock();   
      RecordVtxOccurTimes();       
      BLS();      
      finishing_time = (double)clock(); 
      
      CollectedRev[MultiStartTimes-1]=GetSolutionRevenue();
      UseTime[MultiStartTimes-1]= finishing_time-starting_time;
      
      printf("\n The %d time to run BLS, collected revenue: %d\n", MultiStartTimes, CollectedRev[MultiStartTimes-1]);
      
      CurSltRevenue=GetSolutionRevenue();
      TotalRevenue += CurSltRevenue;
      if(CurSltRevenue > GlobalBestRevenue)
      {
         GlobalBestRevenue=CurSltRevenue;
         SaveCurSltToGlobalBestArray();         
      }        
   }
   
   RestoreSolutionFromGlobalBestArray();  
}//End Multi_BLS()

void RunBLS()
{  
   Multi_BLS();       
   
   FILE *fp; 
   fp = fopen(OutputRltFile, "a+");   
   fprintf(fp, "\n\n %d %d %d %d %d %d %d %d %d %d %d %d %.0f %d", -2, InstCls, InstNo, InstType, Max_Vtx, VtxWithPlusRvnNum, Max_Edge, GetTotalRevenue(), GetTotalBudget()/100, BudgetDivideNum, MaxHop-1, ExpectedRevenue, PrepareTime, RunBLSTimes);
  
   printf("\n\n %d %d %d %d %d %d %d %d %d %d %d %d %.0f %d", -2, InstCls, InstNo, InstType, Max_Vtx, VtxWithPlusRvnNum, Max_Edge, GetTotalRevenue(), GetTotalBudget()/100, BudgetDivideNum, MaxHop-1, ExpectedRevenue, PrepareTime, RunBLSTimes);
  
   for(int i=0;i<RunBLSTimes;i++) 
   {   
     fprintf(fp, " %d %.0f", CollectedRev[i], UseTime[i]); 
     printf(" %d %.0f", CollectedRev[i], UseTime[i]);  
   }   
   fprintf(fp, " %d",-1);  
   fclose(fp);
   
   cout<<endl<<endl;
   return; 
}//End RunBLS()

void SetSumResultFile()
{   
   int len=0;
   OutputRltFile[len++]='A';
   OutputRltFile[len++]='-';
   OutputRltFile[len++]='S';
   OutputRltFile[len++]='u';
   OutputRltFile[len++]='m';
   OutputRltFile[len++]='R';
   OutputRltFile[len++]='e';
   OutputRltFile[len++]='s';
   OutputRltFile[len++]='u';
   OutputRltFile[len++]='l';
   OutputRltFile[len++]='t';
   OutputRltFile[len++]='s';
   OutputRltFile[len++]='-';
   if(Series==1)
     OutputRltFile[len++]='B';
   else
   {
     OutputRltFile[len++]='C';   
      OutputRltFile[len++]='0'+Series-1; 
   }
         
   OutputRltFile[len++]='-';
   OutputRltFile[len++]='V';
   if(Version<0 && Version>=100)
     OutputRltFile[len++]='X';
   else if(Version<10)
     OutputRltFile[len++]='0'+Version;
   else if(Version<100)
   {
     OutputRltFile[len++]='0'+Version/10;  
     OutputRltFile[len++]='0'+Version%10; 
   }
   OutputRltFile[len++]='.';
   OutputRltFile[len++]='t';
   OutputRltFile[len++]='x';
   OutputRltFile[len++]='t';
   OutputRltFile[len++]='\0'; 
}//End ReadParameters()

// Read input parameters of a number of instances
void ReadInstParam()
{
   ifstream FIC;
   FIC.open(InputInstFile);  
  
   if (FIC.fail())
     exit(0);       
      
   FIC >> InstNum;     
   cout << "Number of Instances = " <<InstNum <<endl;     
   
   FILE *fp;      
   fp = fopen(OutputRltFile, "w+");   
   if(!fp)
   {
     cout<<"Fail to save result file to"<<OutputRltFile<<endl;
     getchar();
   }
   fprintf(fp, "%d \n", InstNum);   
   fclose(fp);  
   
   int TempInstCls;
   int TempInstNo;
   int TempInstType;
   int TempBudgetDivideNum;
   int TempHop;     
   for(int i=0;i<InstNum;i++)
   {
      FIC>>TempInstCls>>TempInstNo>>TempInstType>>TempBudgetDivideNum>>TempHop;
      InputParam[i][0]=TempInstCls;
      InputParam[i][1]=TempInstNo;
      InputParam[i][2]=TempInstType;
      InputParam[i][3]=TempBudgetDivideNum;
      InputParam[i][4]=TempHop+1;  
   }    
}// End  ReadInstParam()

void SetInputFileName(int InstIndex)
{
  InstCls=InputParam[InstIndex][0];
  InstNo=InputParam[InstIndex][1];
  InstType=InputParam[InstIndex][2];
      
  int len=0;
  if(InstCls==1) 
  {
    inputfilename[len++]='M';
    inputfilename[len++]='s';
    inputfilename[len++]='t';
    inputfilename[len++]='e';
    inputfilename[len++]='i';
    inputfilename[len++]='n';
    inputfilename[len++]='b';
    if(InstNo<10)  
      inputfilename[len++]='0'+InstNo; 
    else
    {
      inputfilename[len++]='0'+InstNo/10; 
      inputfilename[len++]='0'+InstNo%10;  
    }  
    inputfilename[len++]='.';
    inputfilename[len++]='t';
    inputfilename[len++]='x';
    inputfilename[len++]='t';
    inputfilename[len++]='\0';                  
  }     
  else if(InstCls==2)
  {
    inputfilename[len++]='C';          
    if(InstNo<10)  
      inputfilename[len++]='0'+InstNo; 
    else
    {
      inputfilename[len++]='0'+InstNo/10; 
      inputfilename[len++]='0'+InstNo%10;  
    }      
      inputfilename[len++]='-';
      inputfilename[len++]='1';
      inputfilename[len++]='0';
      if(InstType==100)
        inputfilename[len++]='0';
      
      inputfilename[len++]='.';
      inputfilename[len++]='d';
      inputfilename[len++]='a';
      inputfilename[len++]='t';
      inputfilename[len++]='\0';        
  }
}//End SetInputFileName()

void SetOutfileName(int InstIndex)
{
  BudgetDivideNum=InputParam[InstIndex][3];           
  MaxHop=InputParam[InstIndex][4];      
  cout<<"Input file:"<<inputfilename<<" b="<<BudgetDivideNum<<" h="<<MaxHop-1<<endl;
            
  int len = 0;
  if(InputParam[InstIndex][0]==1)
    outfilename[len++]='B';
  else   
    outfilename[len++]='C';    
        
  if(InputParam[InstIndex][1]<10)  
    outfilename[len++]='0'+InputParam[InstIndex][1]; 
  else
  {
    outfilename[len++]='0'+InputParam[InstIndex][1]/10; 
    outfilename[len++]='0'+InputParam[InstIndex][1]%10;  
  } 
  if(InputParam[InstIndex][2]>=10)  
  {          
    outfilename[len++]='-';
    outfilename[len++]='1';
    outfilename[len++]='0';
    if(InputParam[InstIndex][2]==100)
      outfilename[len++]='0';
  }
      
  outfilename[len++]='(';
  outfilename[len++]='b';
  outfilename[len++]='=';
  if(BudgetDivideNum<10)
    outfilename[len++]='0'+BudgetDivideNum;
  else if(BudgetDivideNum<100)
  {
    outfilename[len++]='0'+BudgetDivideNum/10;
    outfilename[len++]='0'+BudgetDivideNum%10;
  }
  else if(BudgetDivideNum<1000)
  {
    outfilename[len++]='0'+BudgetDivideNum/100;
    outfilename[len++]='0'+(BudgetDivideNum/10)%10;
    outfilename[len++]='0'+BudgetDivideNum%10;
  }
  else if(BudgetDivideNum<10000)
  {
    outfilename[len++]='0'+BudgetDivideNum/1000;
    outfilename[len++]='0'+(BudgetDivideNum/100)%10;
    outfilename[len++]='0'+(BudgetDivideNum/10)%10;
    outfilename[len++]='0'+BudgetDivideNum%10;
  }
  else if(BudgetDivideNum<100000)
  {
    outfilename[len++]='0'+BudgetDivideNum/10000;
    outfilename[len++]='0'+(BudgetDivideNum/1000)%10;
    outfilename[len++]='0'+(BudgetDivideNum/100)%10;
    outfilename[len++]='0'+(BudgetDivideNum/10)%10;
    outfilename[len++]='0'+BudgetDivideNum%10;
  }
  else
  {
    outfilename[len++]='n';
    outfilename[len++]='u';
    outfilename[len++]='l';
    outfilename[len++]='l';
  }
      
  outfilename[len++]=' ';
  outfilename[len++]='h';
  outfilename[len++]='=';
  if(MaxHop-1<10)
    outfilename[len++]='0'+MaxHop-1;
  else if(MaxHop-1<100)
  {
    outfilename[len++]='0'+(MaxHop-1)/10;
    outfilename[len++]='0'+(MaxHop-1)%10;
  }
      
  outfilename[len++]=')';      
  outfilename[len++]='.';
  outfilename[len++]='o';
  outfilename[len++]='u';      
  outfilename[len++]='t';
  outfilename[len++]='\0';     
}// End SetOutfileName()

void OutputSolution()
{
   cout<<"output to file:"<<outfilename<<endl;
   FILE *fp; 
   fp = fopen(outfilename, "w+"); 
   int Edge_Num=0;
   for(int i=0; i<Max_Vtx; i++)
   {
      if(FinalSolution[i]!=Null)
        Edge_Num++;              
   } 
   fprintf(fp, "%d \n", Edge_Num);
   
   for(int i=0; i<Max_Vtx; i++)
   {
      if(FinalSolution[i]!=Null)
        fprintf(fp, "%d %d\n",FinalSolution[i]+1,i+1);              
   }    
         
   fprintf(fp, "Total revenue: %d\n",GetTotalRevenue()); 
   fprintf(fp, "Total edge cost: %d\n",GetTotalBudget()/100); 
   fprintf(fp, "Limited budget: %d\n",Max_Budget/100); 
   fprintf(fp, "Collected revenue: %d\n",GetSolutionRevenue()); 
   fprintf(fp, "Consumed budget: %d\n",GetSolutionBudget()/100); 
    
   fclose(fp) ;
   return;  
}//End OutputSolution()

void RunOneInstance()
{   
   IfDisplayDetail=false;      

   Initializing();       
   Max_Budget = GetTotalBudget()/BudgetDivideNum;  
   GetShortestPathWithHop(MaxHop);
   SaveShortestPath();
   GetCompleteShortestPathWithHop(MaxHop);

   int ReachableRevenue=0;
   for(int i=0;i<Max_Vtx;i++)
   {
      if(CompleteShortestPath[i][0][MaxHop-1]<InfiniteEdgeCost)
        ReachableRevenue+=VtxRevenue[i];                  
   }
   ExpectedRevenue=ReachableRevenue;   
     
   RunBLS(); 
   
   OutputSolution();      
   ReleaseMemory();    
}// End RunOneInstance()

// Run a number of instances sequently
void AutoRunInstancesInBatch()
{
   ReadInstParam();
   for(int i=0;i<InstNum;i++)
   {
      SetInputFileName(i); 
      SetOutfileName(i);        
      // Independently call BLS 10 times to run one instance    
      RunOneInstance();
    } 
}// End AutoRunInstancesInBatch()

int main()
{ 
   SetSumResultFile();    
   AutoRunInstancesInBatch();  
   getchar(); 
   return 0;  
}
