//============================================================================================
// Filename    : ILSMinDiff.cpp
// Author      : Yangming Zhou (yangming@univ-angers.fr)
// Organization: LERIA, Université d’Angers, 2 Boulevard Lavoisier, 49045 Angers, France
// Revised	   : June 2016
// Copyright   : http://dx.doi.org/10.1016/j.knosys.2017.03.028
// Description : An iterated local search algorithm for minimum different dispersion problem
//============================================================================================

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

// Self-definition functions (SDF)
#define SDF_MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
#define SDF_MAX(X,Y)  ((X) < (Y) ? (Y) : (X))
#define SDF_ABS(X)  (((X) < 0) ? -(X) : (X))
#define MAXINT 2147483647
#define MININT -2147483648
#define EPSILON 0.000001
#define NOE 1				// number of elite solutions

char homeDirectory[]="./";	// home directory
char *instanceName; 		// instance file name
char *suffix;
char instancePath[50];		// instance file path
char resultPath[50]; 		// result file path
char statisticPath[50];		// statistical file path

double **diversity;   		// diversity matrix
int nbrTotalVertex;			// total number of elements
int nbrChooseVertex;		// number of chosen elements
int nbrDeltaVertex;
int allocateSpace;
double maxDiversity;
double maxRunTime;			// time limit

double *gain;
int *offspring;
int *bestSol;
int *improvedSol;
int *isChoose;
double bestCost,bestTime;
double improvedCost,improvedTime;
int weakPerturbFactor;
int strongPerturbFactor;
int searchDepth;
double theta;

using namespace std;

/* Read instance */
void read_instance()
{
	ifstream FIC;
    FIC.open(instancePath);
    if(FIC.fail())
    {
    	cout << "### Fail to open file: " << instanceName << endl;
        getchar();
        exit(0);
    }
    if(FIC.eof())
    {
    	cout << "### Fail to open file: " << instanceName << endl;
        exit(0);
    }

    int nbr_pairs = -1;
	FIC >> nbrTotalVertex >> nbrChooseVertex;
	nbrDeltaVertex = nbrTotalVertex-nbrChooseVertex;
	nbr_pairs =(nbrTotalVertex*(nbrTotalVertex-1))/2;

	bestSol = new int[nbrChooseVertex];
	improvedSol = new int[nbrChooseVertex];
	isChoose = new int[nbrTotalVertex];
	gain = new double[nbrTotalVertex];

    diversity = new double*[nbrTotalVertex];
    for(int x = 0; x < nbrTotalVertex; x++)
        diversity[x] = new double[nbrTotalVertex];
    for(int x = 0; x < nbrTotalVertex; x++)
    	for(int y = 0; y < nbrTotalVertex; y++)
    		diversity[x][y] = 0.0;

    double density = 1.0;
    double percent = (double)nbrChooseVertex/(double)(nbrTotalVertex);
    cout << "The statistics of the instance " << instanceName << endl;
    cout << "n = " << nbrTotalVertex << ", m = " << nbrChooseVertex << ", m/n = " << percent << ", and density = " << density << endl;

    maxDiversity = 0.0;
    int x1, x2;
    double x3;
    for(int i = 0; i < nbr_pairs; i++)
    {
        FIC >> x1 >> x2 >> x3;
        if ( x1 < 0 || x2 < 0 || x1 >= nbrTotalVertex || x2 >= nbrTotalVertex )
        {
            cout << "### Read Data Error : line = "<< x1 << ", column = " << x2 << endl;
            exit(0);
        }
    	diversity[x1][x2] = diversity[x2][x1] = x3;
    	if(diversity[x1][x2] > maxDiversity)
    		maxDiversity = diversity[x1][x2];
    }

    cout << "Finish loading data!" << endl;
    FIC.close();
}

/* Verify and store the results */
void check_and_store_result(int *sol,double sol_cost,double sol_time)
{
	FILE *out;
	double max_gain = -1.0;
	double min_gain = (double)MAXINT;
	double sol_true_cost;

	for(int i = 0; i < nbrChooseVertex; i++)
	{
		gain[sol[i]] = 0.0;
		for(int j = 0; j < nbrChooseVertex; j++)
			if(sol[j] != sol[i])
				gain[sol[i]] += diversity[sol[i]][sol[j]];

		if(gain[sol[i]] > max_gain)
			max_gain = gain[sol[i]];
		if(gain[sol[i]] < min_gain)
			min_gain = gain[sol[i]];
	}
	sol_true_cost = max_gain - min_gain;

	if(SDF_ABS(sol_true_cost-sol_cost) > EPSILON)
	{
		printf("Find a error solution!\n");
		printf("sol_cost = %lf, while sol_true_cost = %lf\n",sol_cost,sol_true_cost);
		exit(0);
	}

	out = fopen(resultPath, "a+");
	// Printing parameters, data characteristics and some statistical results
	fprintf(out,"   Statistical information:                                 \n");
	fprintf(out,"   Total number of elements (n)               	= %d\n",nbrTotalVertex);
	fprintf(out,"   Number of selected elements (m)           	= %d\n",nbrChooseVertex);
	fprintf(out,"   Weak perturbation factor                	= %d\n",weakPerturbFactor);
	fprintf(out,"   Strong perturbation factor                	= %d\n",strongPerturbFactor);
	fprintf(out,"   theta                						= %lf\n",theta);
	fprintf(out,"   Search depth                				= %d\n",searchDepth);
	fprintf(out,"   Maximum diversity                			= %lf\n",maxDiversity);
	fprintf(out,"   Max running time (t_max,second)             = %lf\n",maxRunTime);
	fprintf(out,"   Found best solution:\n");
	for(int i = 0; i < nbrChooseVertex; i++)
		fprintf(out, "%d ", sol[i]);
	fprintf(out, "\n");
	fprintf(out,"   best solution = %lf, best time = %lf\n",sol_cost,sol_time);
	fclose(out);
}

/* local search with best improvement */
void descent_based_search(int *sol,double max_time,double begin_time)
{
	int *sol1;
	int *add_elements;
	int *address;
	double est_gain;
	double est_to_add_gain;
	double est_max_gain,est_min_gain;
	double delta_cost;
	int to_remove_neighbor[400];
	int to_add_neighbor[400];
	int nbr_best_neighbor;
	double best_delta_cost;
	int index;
	int to_remove,to_add;
	int ad0,ad1;
	int is_improve = 1;
	double sol_cost;
	double max_gain,min_gain;
	sol1 = new int[nbrChooseVertex];
	add_elements = new int[nbrDeltaVertex];
	address = new int[nbrTotalVertex];

	for(int i = 0; i < nbrChooseVertex; i++)
		improvedSol[i] = sol[i];
	improvedTime = (clock()-begin_time)/CLOCKS_PER_SEC;

	min_gain = (double)MAXINT;
	max_gain = -1.0;
	for(int i = 0; i < nbrChooseVertex; i++)
	{
		gain[sol[i]] = 0.0;
		for(int j = 0; j < nbrChooseVertex; j++)
			if(sol[j] != sol[i])
				gain[sol[i]] += diversity[sol[i]][sol[j]];

		if(gain[sol[i]] > max_gain)
			max_gain = gain[sol[i]];
		if(gain[sol[i]] < min_gain)
			min_gain = gain[sol[i]];
	}
	sol_cost = max_gain - min_gain;
	improvedCost = sol_cost;

	memset(isChoose,0,allocateSpace);
	memset(address,0,allocateSpace);
	for(int i = 0; i < nbrChooseVertex; i++)
		isChoose[sol[i]] = 1;

	int x1,x2;
	x1 = x2 = 0;
	for(int i = 0; i < nbrTotalVertex; i++)
	{
		if(isChoose[i] == 1)
		{
			sol1[x1] = i;
			address[i] = x1;
			x1++;
		}
		else if(isChoose[i] == 0)
		{
			add_elements[x2] = i;
			address[i] = x2;
			x2++;
		}
		else
		{
			cout << "input solution is not correct: " << "is_choose[" << i << "] =" << isChoose[i] << endl;
			exit(0);
		}
	}

	// Search until a local minimum is reached
	while(is_improve != 0)
	{
		is_improve = 0;
		best_delta_cost = (double)MAXINT;
		nbr_best_neighbor = 0;
		// find the best neighboring moves
		for(int i = 0; i < nbrChooseVertex; i++)
			for(int j = 0; j < nbrDeltaVertex; j++)
			{
				to_remove = sol1[i];
				to_add = add_elements[j];

				// calculate the move gain
				est_max_gain = -1.0;
				est_min_gain = (double)MAXINT;
				est_to_add_gain = 0.0;
				for(int k = 0; k < nbrChooseVertex; k++)
					if(sol1[k] != to_remove)
					{
						est_gain = gain[sol1[k]] - diversity[sol1[k]][to_remove] + diversity[sol1[k]][to_add];
						if(est_gain > est_max_gain)
							est_max_gain = est_gain;
						if(est_gain < est_min_gain)
							est_min_gain = est_gain;

						// estimation gain of adding element v
						est_to_add_gain += diversity[to_add][sol1[k]];
					}

				// move gain;
				delta_cost = SDF_MAX(est_to_add_gain,est_max_gain)-SDF_MIN(est_to_add_gain,est_min_gain)-sol_cost;

				if(delta_cost < best_delta_cost)
				{
					best_delta_cost = delta_cost;
					to_remove_neighbor[0] = to_remove;
					to_add_neighbor[0] = to_add;
					nbr_best_neighbor = 1;
				}
				else if((SDF_ABS(delta_cost-best_delta_cost) < EPSILON) && (nbr_best_neighbor < 400))
				{
					to_remove_neighbor[nbr_best_neighbor] = to_remove;
					to_add_neighbor[nbr_best_neighbor] = to_add;
					nbr_best_neighbor++;
				}
			}

		// select a best move and perform a move
		if(best_delta_cost < 0.0 && nbr_best_neighbor > 0)
		{
			// select a move
			index = rand()%nbr_best_neighbor;
			to_remove = to_remove_neighbor[index];
			to_add = to_add_neighbor[index];

			is_improve = 1;
			sol_cost += best_delta_cost;

			// make a move
	        ad0 = address[to_remove];
	        ad1 = address[to_add];
	        sol1[ad0] = to_add;
	        address[to_add] = ad0;
	        add_elements[ad1] = to_remove;
	        address[to_remove] = ad1;

			// update after making a move
			gain[to_add] = 0.0;
			for(int i = 0; i < nbrChooseVertex; i++)
				if(sol1[i] != to_add)
				{
					gain[sol1[i]] += -diversity[sol1[i]][to_remove] + diversity[sol1[i]][to_add];
					gain[to_add] += diversity[to_add][sol1[i]];
				}
		}

		if((clock()-begin_time)/CLOCKS_PER_SEC > max_time)
			break;
	}

	improvedTime = (clock()-begin_time)/CLOCKS_PER_SEC;
	improvedCost = sol_cost;
	for(int i = 0; i < nbrChooseVertex; i++)
		improvedSol[i] = sol1[i];

	delete [] sol1;
	delete [] add_elements;
	delete [] address;
}

void generate_an_elite_solution(double max_time,double begin_time)
{
    int *add_elements;
    int *initial_sol;
    int nbr_to_add_element;
    int nbr_element;
    int nbr_sol;
    int index;
    initial_sol = new int[nbrChooseVertex];
    add_elements = new int[nbrTotalVertex];

    bestCost = (double)MAXINT;
    nbr_sol = 0;
    while(nbr_sol < NOE)
    {
		nbr_to_add_element = 0;
		for(int j = 0; j < nbrTotalVertex; j++)
			add_elements[nbr_to_add_element++] = j;

		memset(isChoose,0,allocateSpace);
		nbr_element = 0;
		while(nbr_element < nbrChooseVertex)
		{
			index = rand()%nbr_to_add_element;
			initial_sol[nbr_element] = add_elements[index];
			isChoose[add_elements[index]] = 1;
			nbr_element++;

			nbr_to_add_element--;
			add_elements[index] = add_elements[nbr_to_add_element];
		}

        descent_based_search(initial_sol,max_time,begin_time);
        //tabu_search(initial_sol,max_time,begin_time);

        if(improvedCost < bestCost)
        {
        	bestCost = improvedCost;
        	bestTime = improvedTime;
        	for(int i = 0; i < nbrChooseVertex; i++)
        		bestSol[i] = improvedSol[i];
        }
        nbr_sol++;
   }

    delete [] initial_sol;
    delete [] add_elements;
}

void weak_perturb_operation(int *sol,int *best_sol,double best_sol_cost)
{
	int *add_elements;
	int *address;
	int nbr_perturb;
	int nbr_pick_time;
	int index_add1,index_remove1;
	int to_add1,to_remove1;
	int index_add2,index_remove2;
	int to_add2,to_remove2;
	double delta_cost1,delta_cost2;
	double sol_cost;
	double est_max_gain,est_min_gain;
	double est_gain,est_to_add_gain;
	int ad0,ad1;
	add_elements = new int[nbrDeltaVertex];
	address = new int[nbrTotalVertex];

    for(int i = 0; i < nbrChooseVertex; i++)
    {
    	gain[best_sol[i]] = 0.0;
        for(int j = 0; j < nbrChooseVertex; j++)
        	if(best_sol[j] != best_sol[i])
        		gain[best_sol[i]] += diversity[best_sol[i]][best_sol[j]];
    }
    sol_cost = best_sol_cost;

	memset(isChoose,0,allocateSpace);
	for(int i = 0; i < nbrChooseVertex; i++)
		isChoose[best_sol[i]]++;

	int x1,x2;
	x1 = x2 = 0;
	for(int i = 0; i < nbrTotalVertex; i++)
	{
		if(isChoose[i] == 1)
		{
			sol[x1] = i;
			address[i] = x1;
			x1++;
		}
		else if(isChoose[i] == 0)
		{
			add_elements[x2] = i;
			address[i] = x2;
			x2++;
		}
		else
		{
			cout << "input solution is not correct: " << "is_choose[" << i << "] =" << isChoose[i] << endl;
			exit(0);
		}
	}

	nbr_perturb = 0;
	while(nbr_perturb < weakPerturbFactor)
	{
		// randomly select a neighboring solution
		index_add1 = rand()%nbrDeltaVertex;
		index_remove1 = rand()%nbrChooseVertex;
		to_add1 = add_elements[index_add1];
		to_remove1 = sol[index_remove1];

		// calculate delta cost
		est_max_gain = -1.0;
		est_min_gain = (double)MAXINT;
		est_to_add_gain = 0.0;
		for(int k = 0; k < nbrChooseVertex; k++)
			if(sol[k] != to_remove1)
			{
				est_gain = gain[sol[k]] - diversity[sol[k]][to_remove1] + diversity[sol[k]][to_add1];

				if(est_gain > est_max_gain)
					est_max_gain = est_gain;
				if(est_gain < est_min_gain)
					est_min_gain = est_gain;

				// estimation gain of adding element v
				est_to_add_gain += diversity[sol[k]][to_add1];
			}
		delta_cost1 = SDF_MAX(est_to_add_gain,est_max_gain)-SDF_MIN(est_to_add_gain,est_min_gain)-sol_cost;

		nbr_pick_time = 0;
		while(nbr_pick_time < nbrTotalVertex)
		{
			// randomly pick a neighboring solution
			index_add2 = rand()%nbrDeltaVertex;
			index_remove2 = rand()%nbrChooseVertex;
			to_add2 = add_elements[index_add2];
			to_remove2 = sol[index_remove2];

			// calculate delta cost
			est_max_gain = -1.0;
			est_min_gain = (double)MAXINT;
			est_to_add_gain = 0.0;
			for(int k = 0; k < nbrChooseVertex; k++)
				if(sol[k] != to_remove1)
				{
					est_gain = gain[sol[k]] - diversity[sol[k]][to_remove1] + diversity[sol[k]][to_add1];

					if(est_gain > est_max_gain)
						est_max_gain = est_gain;
					if(est_gain < est_min_gain)
						est_min_gain = est_gain;

					// estimation gain of adding element v
					est_to_add_gain += diversity[sol[k]][to_add1];
				}
			delta_cost2 = SDF_MAX(est_to_add_gain,est_max_gain)-SDF_MIN(est_to_add_gain,est_min_gain)-sol_cost;

			// record the best move
			if(delta_cost2 < delta_cost1)
			{
				delta_cost1 = delta_cost2;
				index_add1 = index_add2;
				index_remove1 = index_remove2;
				to_add1 = to_add2;
				to_remove1 = to_remove2;
			}
			nbr_pick_time++;
		}

		// make the move
        ad0 = address[to_remove1];
        ad1 = address[to_add1];
        sol[ad0] = to_add1;
        address[to_add1] = ad0;
        add_elements[ad1] = to_remove1;
        address[to_remove1] = ad1;

        // update the gain;
		gain[to_add1] = 0.0;
		for(int i = 0; i < nbrChooseVertex; i++)
			if(sol[i] != to_add1)
			{
				gain[sol[i]] += -diversity[sol[i]][to_remove1] + diversity[sol[i]][to_add1];
				gain[to_add1] += diversity[to_add1][sol[i]];
			}

		nbr_perturb++;
	}

	delete [] add_elements;
	delete [] address;
}

void local_optima_exploring(int *sol,int search_depth,double max_time,double begin_time)
{
	int count;
	int *best_sol;
	double best_sol_cost;
	double best_sol_time;
	best_sol = new int[nbrChooseVertex];

	/* local search */
	descent_based_search(sol,max_time,begin_time);

	best_sol_cost = improvedCost;
	best_sol_time = improvedTime;
	for(int i = 0; i < nbrChooseVertex; i++)
		best_sol[i] = improvedSol[i];

	count = 0;
	while(count < search_depth)
	{
		/* weak perturb */
		weak_perturb_operation(sol,best_sol,best_sol_cost);
		/* local search */
		descent_based_search(sol,max_time,begin_time);

		if(improvedCost < best_sol_cost)
		{
			best_sol_cost = improvedCost;
			best_sol_time = improvedTime;
			for(int i = 0; i < nbrChooseVertex; i++)
				best_sol[i] = improvedSol[i];
			count = 0;
		}
		else
			count++;
	}

	// return best solution
	improvedCost = best_sol_cost;
	improvedTime = best_sol_time;
	for(int i = 0; i < nbrChooseVertex; i++)
		improvedSol[i] = best_sol[i];

	delete [] best_sol;
}

void local_optima_escaping(int *sol,int perturb_strength)
{
	int *add_elements;
	int *address;
	int nbr_perturb;
	int index_add,index_remove;
	int to_add,to_remove;
	int ad0,ad1;
	add_elements = new int[nbrDeltaVertex];
	address = new int[nbrTotalVertex];

	memset(isChoose,0,allocateSpace);
	memset(address,0,allocateSpace);
	for(int i = 0; i < nbrChooseVertex; i++)
		isChoose[sol[i]]++;

	int x1,x2;
	x1 = x2 = 0;
	for(int i = 0; i < nbrTotalVertex; i++)
	{
		if(isChoose[i] == 1)
		{
			sol[x1] = i;
			address[i] = x1;
			x1++;
		}
		else if(isChoose[i] == 0)
		{
			add_elements[x2] = i;
			address[i] = x2;
			x2++;
		}
		else
		{
			cout << "input solution is not correct: " << "is_choose[" << i << "] =" << isChoose[i] << endl;
			exit(0);
		}
	}

	nbr_perturb = 0;
	while(nbr_perturb < perturb_strength)
	{
		index_add = rand()%nbrDeltaVertex;
		to_add = add_elements[index_add];
		index_remove = rand()%nbrChooseVertex;
		to_remove = sol[index_remove];

		// make the move
        ad0 = address[to_remove];
        ad1 = address[to_add];
        sol[ad0] = to_add;
        address[to_add] = ad0;
        add_elements[ad1] = to_remove;
        address[to_remove] = ad1;

		nbr_perturb++;
	}

	delete [] add_elements;
	delete [] address;
}

void ILSMinDiff(char *result_file,double max_time)
{
	int iter;
	int *sol;
	int no_improve_iter;
	double begin_time = clock();
	double theta;
	sol = new int[nbrChooseVertex];

	// search depth of local optima exploring phase
	searchDepth = 5;
	//weak perturbation factor
	if(nbrTotalVertex < 500 || (nbrTotalVertex == 500 && nbrTotalVertex/nbrChooseVertex < 10))
		weakPerturbFactor = 3;
	else
	weakPerturbFactor = 2;

	// strong perturbation factor
	theta = 1.0;
	strongPerturbFactor = (int)(theta*(nbrTotalVertex/nbrChooseVertex));

	/* generate a good initial solution */
	generate_an_elite_solution(max_time,begin_time);

	for(int i = 0; i < nbrChooseVertex; i++)
		sol[i] = bestSol[i];

	cout << "used time = " << (clock()-begin_time)/CLOCKS_PER_SEC << endl;
	/* iteratively run ITS */
	iter = 0;
	printf("iter = %d,best cost = %.3lf, and improved cost = %.3lf\n",iter,bestCost,improvedCost);
	while(((clock()-begin_time)/CLOCKS_PER_SEC) < max_time)
	{
		local_optima_exploring(sol,searchDepth,max_time,begin_time);

        //record the best solution
        if(improvedCost < bestCost)
        {
        	bestTime = improvedTime;
        	bestCost = improvedCost;
        	for(int i = 0; i < nbrChooseVertex; i++)
        		bestSol[i] = improvedSol[i];
        	no_improve_iter = 0;
        }
        else
        	no_improve_iter++;

        local_optima_escaping(improvedSol,strongPerturbFactor);

        for(int i = 0; i < nbrChooseVertex; i++)
    		sol[i] = improvedSol[i];

        iter++;
        printf("iter = %d,best cost = %.3lf, and improved cost = %.3lf\n",iter,bestCost,improvedCost);
	}

	delete [] sol;
}

/* Main procedure */
int main(int argc,char *argv[])
{
	FILE *statistic;
	int nbr_repeat;
	int count;
	char numbers[41][2]={{'0'},{'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'},{'8'},{'9'},{'1','0'},
			{'1','1'},{'1','2'},{'1','3'},{'1','4'},{'1','5'},{'1','6'},{'1','7'},{'1','8'},{'1','9'},{'2','0'},
			{'2','1'},{'2','2'},{'2','3'},{'2','4'},{'2','5'},{'2','6'},{'2','7'},{'2','8'},{'2','9'},{'3','0'},
			{'3','1'},{'3','2'},{'3','3'},{'3','4'},{'3','5'},{'3','6'},{'3','7'},{'3','8'},{'3','9'},{'4','0'}};
	double sd;
	double avg_sol_cost;
	double avg_sol_time;
	double global_best_sol_cost;
	double global_worst_sol_cost;
	int nbr_best_sol;

	if(argc == 3)
    {
		instanceName = argv[1]; // instance name
		nbr_repeat = atoi(argv[2]); // total number of trials
    }
    else
    {
        cout << endl << "### Input the following parameters ###" << endl;
        cout << "Instance file name, and number of repeats" << endl;
        exit(0);
    }

    // Instance and statistical file path
    strcpy(instancePath, homeDirectory);
    strcat(instancePath, "instances/MDG-c/");
    strcat(instancePath, instanceName);

    strcpy(statisticPath, homeDirectory);
    strcat(statisticPath, "statistics/MDG-c/");
    strcat(statisticPath, instanceName);
    strcat(statisticPath, ".statistic");

	if((statistic = fopen(statisticPath,"w")) == NULL)
    {
		printf("Open failed for output  %s\n",statisticPath);
		exit(1);
    }

    //Read the instance
    read_instance();

    double best_sol_time[nbr_repeat];
    double best_sol_cost[nbr_repeat];
    allocateSpace = nbrTotalVertex*sizeof(int);

	// Repeat multiple runs
	maxRunTime = nbrTotalVertex;
	fprintf(statistic,"Maximum allowable running time = %.3lf\n",maxRunTime);
	fprintf(statistic,"---------------------------------------------------------\n");
	for(int i = 1; i <= nbr_repeat; i++)
    {
		srand((unsigned) time(NULL));
	    strcpy(resultPath, homeDirectory); // result file path
	    strcat(resultPath, "results/MDG-c/");
	    strcat(resultPath, instanceName);
	    strcat(resultPath, ".result");

		count = strlen(resultPath);
		if(i < 10)
		{
			resultPath[count] = numbers[i][0];
		}
		else
		{
			resultPath[count] = numbers[i][0];
			count++;
			resultPath[count] = numbers[i][1];
		}
		count++;
		resultPath[count]='\0';

		// Run the ILSMinDiff algorithm
		ILSMinDiff(resultPath,maxRunTime);

		// Check and store the results
		check_and_store_result(bestSol,bestCost,bestTime);

		// Statistical results
		best_sol_cost[i-1] = bestCost;
		best_sol_time[i-1] = bestTime;
		fprintf(statistic,"%.3lf    %.3lf    %d\n",bestCost,bestTime,i);

    }

	// Compute the statistical results
	global_best_sol_cost = (double)MAXINT;
	global_worst_sol_cost = -1.0;
	avg_sol_cost = 0.0;
	avg_sol_time = 0.0;
	for(int i = 0; i < nbr_repeat; i++)
	{
		avg_sol_cost += best_sol_cost[i];
		avg_sol_time += best_sol_time[i];
		// find out the best solution
		if(best_sol_cost[i] < global_best_sol_cost + EPSILON)
		{
			global_best_sol_cost = best_sol_cost[i];
			nbr_best_sol = 1;
		}
		else if(SDF_ABS(best_sol_cost[i]-global_best_sol_cost) < EPSILON)
			nbr_best_sol++;
		// find out the worst solution
		if(best_sol_cost[i] > global_worst_sol_cost + EPSILON)
			global_worst_sol_cost = best_sol_cost[i];
	}

	avg_sol_cost = avg_sol_cost/(double)nbr_repeat;
	avg_sol_time = avg_sol_time/(double)nbr_repeat;

	// standard deviation
	sd = 0.0;
	for(int i = 0; i < nbr_repeat; i++)
		sd += (best_sol_cost[i] - avg_sol_cost)*(best_sol_cost[i] - avg_sol_cost);
	sd = sd/(double)nbr_repeat;
	if(sd > 0.0)
		sd = sqrt(sd);

	fprintf(statistic,"---------------------------------------------------------\n");
	fprintf(statistic,"%.3lf	%.3lf    %.3lf    %.3lf    %.3lf\n",global_best_sol_cost,avg_sol_cost,global_worst_sol_cost,avg_sol_time,sd);
	cout << "Instance:" << instanceName << ", the best cost = " << global_best_sol_cost << ", average cost = "<< avg_sol_cost << ", and the worst cost = " << global_worst_sol_cost << endl;

    // Free memory
    delete [] gain;
    delete [] isChoose;
    delete [] offspring;
    delete [] bestSol;
    delete [] improvedSol;

    for(int i = 0; i < nbrTotalVertex; i++)
    	delete diversity[i];

    fclose(statistic);
    return 0;

}
