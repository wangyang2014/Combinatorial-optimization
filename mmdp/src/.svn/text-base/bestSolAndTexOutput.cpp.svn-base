#include "bestSolAndTexOutput.h"

//This routine returns the number of iterations with no new better solution found
int itersNoNewBetter(){
    return it-lastIterWithNewBestPRecord;
}

//Return value: 0 -  it is not better, 1 - better, 2=better that all time best
//This routines should work ok if there are less than k vertices acctivated
int checkNewBetterSolution(){

    //Internal technical note:
    //Unless epsilon is used, it can consider that
    //bestMinDist=minDist =>
    //minDist>bestMindist
    int firstCriterionImproved = (minDist>bestMinDist+epsilon);
    if (firstCriterionImproved){
        if(itersNoNewBetter()>maxItersNoNewBetter)
            maxItersNoNewBetter = itersNoNewBetter();
        lastIterWithNewBestPRecord = it;
        timeAtNewBetter = timeCpuPassedMSec()/1000;

    }
    if(firstCriterionImproved||(minDist==bestMinDist&&sumDist>bestSumDist+epsilon)){



        if(bestX==NULL)
            bestX = new int[n];
        bestMinDist = minDist;
        bestSumDist = sumDist;
        for(int ii=0;ii<n;ii++)
            bestX[ii]=x[ii];
        if(minDist>=bestKnown+0.01&&h==k){       //previous records show two decimals, but h==k because we take
            saveBestSolutionBetterPerf();        //into consideration only complete solutions here
            return 2;
        }
        return 1;
    }
    return 0;
}


//Useful when going from h to h+1.
void resetBestSolution(){
    if(bestX==NULL){
        delete bestX;
        bestX=NULL;
    }
    bestMinDist=-1;
    lastIterWithNewBestPRecord=it;
}

//The verification will be done on vector x
//passed as argument to this function
int checkMMDP(int* x,double minDist){
    double minDistCheck=INT_MAX;
    for(int ii=0;ii<n;ii++)
        if(x[ii])
            for(int jj=ii+1;jj<n;jj++)
                if(x[jj])
                    if(m[ii][jj]<minDistCheck)
                        minDistCheck=m[ii][jj];
    //saveSolution();
    return (minDist==minDistCheck);
}
//The verification will be done on array x and double sumDist
//passed as argument to this function
int checkMDP(int*x, double &sumDist){
    double sumDistCheck=0;
    for(int ii=0;ii<n;ii++)
        if(x[ii])
            for(int jj=ii+1;jj<n;jj++)
                if(x[jj])
                    sumDistCheck+=m[ii][jj];
    //cout<<sumDistCheck <<endl;
    //cout<<sumDist<<endl;
    //cout<<(sumDist==sumDistCheck)<<endl;
    //if(sumDist!=sumDistCheck)
    //    cout<<"Warning: sumDist="<<sumDist<<" but sumDistCheck="<<sumDistCheck<<endl;

    double sumEpsilon=0.01;//In long run, rounding errors can accumulate here
                          //
    if(fabs(sumDist-sumDistCheck)<sumEpsilon){
        return 1;
    }
    else
        if(fabs(sumDist-sumDistCheck)<sumEpsilon*10){//Rounding errors should only be corrected
            sumDist=sumDistCheck;
            //You shoud recompute everything
            return 1;
        }
    cerr<<"SumDist diff:"<<sumDist-sumDistCheck<<endl;
    return 0;
    //return (sumDsist==sumDistCheck);

}
void saveThisSpecificSolution(int* x, double minDist, int keepOutputFolderInFilename){
    cout<<"Saving the instance";
    string folderToPutFile="./";

    //Take abc/dsd/dasd from abc/dsd/dasd/
    string trimedOutputFolder = outputFolder.substr(0,outputFolder.find_last_not_of("/ ")+1);

    //currentTopFolder should be dasd
    string currentTopFolder=trimedOutputFolder;
    string::size_type lastSlash = trimedOutputFolder.find_last_of("/");
    if(lastSlash!=string::npos){
        currentTopFolder = trimedOutputFolder.substr(lastSlash+1);
        folderToPutFile = trimedOutputFolder.substr(0,lastSlash)+"/";
    }
    //Don't put the folder name in the filename to be written
    if(!keepOutputFolderInFilename)
        currentTopFolder = "";
    ofstream solOut((folderToPutFile+currentTopFolder+instance+".sol.txt").c_str());
    cout<<" to "<<(folderToPutFile+currentTopFolder+instance+".sol.txt")<<endl;
    if(solOut.is_open()){
        solOut<<"Instance name="<<instance<<endl;
        solOut<<"n="<<n<<endl;
        solOut<<"m="<<k<<endl;
        solOut<<"Previous bestKnown (according to our records)="<<bestKnown<<endl;
        solOut<<"Current solution's selected vertices (0-indexed):"<<endl;
        for(int ii=0;ii<n;ii++)
            if(x[ii])
                solOut<<setw(4)<<ii;
        solOut<<endl;
        solOut<<"This solution has a MMDP value of "<<minDist<<endl;
        printMatrix(solOut,x);
        solOut.close();
    }
    else
        cout<<"Warning: I could not write the solution to file:"<<folderToPutFile+currentTopFolder+instance+".txt"<<endl;

}

void saveThisSpecificSolution(int* x, double minDist){
    saveThisSpecificSolution(x, minDist, 1);
}

void saveBestSolutionBetterPerf(){
    saveThisSpecificSolution(bestX,bestMinDist,1);
}

void saveBestSolutionEqualPerf(){
    saveThisSpecificSolution(bestX,bestMinDist,0);
}


void finishCheckAndWriteTex(int *x, double minDist, double sumDist){
    texFile<<setw(15)<<instance;
    texFile<<setw(1)<<"&";
    texFile<<setw(10)<<n;texFile<<setw(1)<<"&";
    texFile<<setw(10)<<k;texFile<<setw(1)<<"&";
    if(!checkMMDP(x,minDist))
        texFile<<setw(10)<<"MMDP CHECK UNUSUAL"<<endl;
    if(minDist>bestKnown+0.0045){//10.005 could be better than 10.00 (unless 10.00 is rounded from 10.0049999?)
        texFile<<setw(18)<<"\\underline{\\textbf{"<<setw(6)<<keepFourDecimals(minDist)<<setw(2)<<"}}";
        //saveSolution();
    }
    else if(keepTwoDecimals(minDist)==keepTwoDecimals( bestKnown))
        texFile<<setw(8)<<"\\textbf{"<<setw(17)<<keepTwoDecimals(minDist)<<setw(1)<<"}";
    else if(keepTwoDecimals(minDist)==keepTwoDecimals( bestKnown)-0.01)
        texFile<<setw(8)<<"\\textbf{"<<setw(17)<<keepFourDecimals(minDist)<<setw(1)<<"}";
    else if(minDist>bestKnown-0.0055)//if bestknown=10.01 rounded from 10.005, than minDist=10.0049 may be the same
        texFile<<setw(8)<<"\\textbf{"<<setw(17)<<keepFourDecimals(minDist)<<setw(1)<<"}";
    else
        texFile<<setw(26)<<keepTwoDecimals(minDist);
    texFile<<setw(1)<<"&";
    if(!checkMDP(x,sumDist))
        texFile<<setw(10)<<"MDP CHECK UNUSUAL"<<endl;
    texFile<<setw(10)<<keepTwoDecimals(sumDist)<<setw(1)<<"&";
    if(timeCpuPassedMSec()>=1000)
        texFile<<setw(10)<<(timeCpuPassedMSec()/1000);
    else
        texFile<<setw(10)<<"$<1$";
    texFile<<setw(1)<<"&";
    texFile<<setw(10);
    if(bestKnown!=INT_MAX)
        texFile<<bestKnown;
    else
        texFile<<"not available";
    texFile<<setw(1)<<"\\\\";
    texFile<<"%&"<<10000*(keepFourDecimals(1.0-(keepFourDecimals(minDist)/bestKnown)))<<"&";
    //cerr<<(keepTwoDecimals(1.0-((minDist)/bestKnown)))<<endl;
    if(minDist>bestKnown-0.0055)
        texFile<<"1";
    else
        texFile<<"0";
    texFile<<" &"<<maxItersNoNewBetter/n<<"&"<<timeAtNewBetter;
    texFile<<"\n";
    texFile.close();
}

void finishCheckAndWriteTex(){
    finishCheckAndWriteTex(bestX,bestMinDist,bestSumDist);
}
void printBestSolution(){
    cout<<"*****************"<<endl;
    cout<<"Final best="<<bestMinDist<<endl;
}

double bestMinDistUntilNow(){
    return bestMinDist;
}

