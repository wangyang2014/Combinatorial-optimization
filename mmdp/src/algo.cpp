#include "algo.h"
#include "choices.h"

//Above object is useful if ones wants to keep all selected vertices in an ordered list (see paragraph above Section 2.5 in the paper)
//This way, one can go through all selected vertices in O(k) instead of O(n)
//Adding a vertex is done in O(n/k) --- in our implementation (in average)
//Dropping a vertex is done in O(1).
//The wasted O(n/k) of the add operation is compensated by several gains (each iteration) of order (O(n)-O(k))
//We could have use std::set, but it looks quite inefficient to go through these point using (more complex) iterators
//see also "choices.h"
#ifndef DISABLE_FAST_ORDERED_LIST
//set<int>    selected;
#include "fastSet.cpp"
fastSet     selected;
#endif

/*
 * ==================================================================
 *         TABU LIST ROUTINES.
 *         The Tabu lists routines (initially using other parameter values) are now "ultra-simplified" and use the values of
 *               -  0 (not active for addition) or
 *               -  h-1 (selecting the "oldest" element for dropping)
 * ==================================================================
 */

//The variable h below is the number of vertices to be selected. It goes from 1 to k in the constructive phase. It is k in the local search
void setGenericTabuTenure(){
    genericAddTenure  = 0;     //NOW 0 //;min(10,(int)(0.5*(n-NoSelectedVtx)));
    genericDropTenure = h-1;   //NOW h-1, only the last one at level h is OK //min(20,(int)(0.25*((double)NoSelectedVtx)));
}
void computeDropTenure(){
    dropTenure    = genericDropTenure;
//  int deviation = (int)(0.25*((double)dropTenure))*0;
//  dropTenure    = dropTenure+randomInt(-deviation,deviation);
}
void computeAddTenure(){
    //setting the Tabu tenure
    addTenure     = genericAddTenure ;//now 0
    int deviation = 0;                //previously used in some tests, now deviation=0
    //Reactive try 1; using the number of iteration that do not improve upon the best sol. found so far
    //if(it%PRINTRATE==0)cout<<itersNoNewBetter()<<"->"<<0.1*((double)addTenure*itersNoNewBetter()/(n*maxIter) )<<endl;
    //addTenure     = addTenure + 0.1*((double)addTenure*itersNoNewBetter()/(n*maxIter) );
    //Reactive try 2; using the number of cycles with a constant best level during all cycle iterations
    //addTenure    += constantCycles;    //Adding some deviations to this generic Tabu tenure

    deviation = max(deviation,2);       //a val of deviation = 1 would be too small, keeping this instr for exceptions
    addTenure = addTenure+randomInt(-deviation,deviation);
    addTenure = min(addTenure,n-k-1);   //addTenure must not exceed n-k-1 because this would make all moves Tabu

}


/*
 * ==================================================================
 *         Routines to make calculations after adding/dropping some chosen elements
 *            - the parameter vtx is the chosen element to be added/dropped
 * ==================================================================
 */

//Function add and drop change iteration counter and Tabu status;
void add(int vtx){
    assert(!x[vtx]);
    x[vtx]=1;
    #ifndef DISABLE_FAST_ORDERED_LIST
    selected.insert(vtx);
    #endif
    iter[vtx]=it;
    computeDropTenure();
    //dropTenure=k-1;
    dropTabu[vtx]=it+dropTenure;
    addTabu[vtx]=INT_MAX;        // an added point should not be added again
    it++;
    #ifndef NDEBUG
    if(it==0){
        cerr<<"Attention: The number of iteration is higher than INT_MAX, you can try using long data type for it, addTabu, etc.";
        for(int ii=0;ii<n;ii++){
            addTabu[ii]=-1;
            dropTabu[ii]=-1;
        }

    }
    #endif

}
void drop(int vtx){
    //saveThisSpecificSolution(x,minDist);
    assert(x[vtx]);
    x[vtx]=0;
    #ifndef DISABLE_FAST_ORDERED_LIST
    selected.erase(vtx);
    #endif
    iter[vtx]=it;
    computeAddTenure();

    addTabu[vtx]=it+addTenure;
    //cerr<<addTenure<<",";
    //If already dropped, it can not be redropped
    dropTabu[vtx]=INT_MAX;
}



//ROUTINE BELOW CORESPONDS TO case (3) of Algorithm 4 of the paper (A simple and effective algorithm for the MaxMin Diversity Problem)
//This operation takes O(n) with DISABLE_FAST_ORDERED_LIST or O(n/k) otherwise
//This routine is applied when point vtx looses its values minx and nmin.
//Usually because the vertex A that linked to vtx via a minimum distance becomes unsellected;
void recomputeMinDistAndNMin(int vtx){
    recomputed++;
    minx[vtx]=INT_MAX;

    #ifdef DISABLE_FAST_ORDERED_LIST
    for(int ii=0;ii<n;ii++)
        if(x[ii])
    #else
    //Version using sts::set iterators is quite heavy
    //set<int>::iterator xIterator=selected.begin();for(int ii=*xIterator;xIterator!=selected.end();++xIterator,ii=*xIterator)
    for(int ii=selected.firstElem();ii<n;ii=selected.nextElem(ii))
    #endif
            if(ii!=vtx){
                if(m[ii][vtx]<minx[vtx]){
                    minx[vtx]=m[ii][vtx];
                    nmin[vtx]=1;
                }
                else
                if(m[ii][vtx]==minx[vtx])
                    nmin[vtx]++;
            }
}

//See Algorithm 4 for details on Algorithm below
void updateWhenDropped(int oldVtx){
    //The values sumx,minx, sumDist for all vertices could be affected after dropping oldVtx
    //Except oldVtx is unchanged
    sumDist-=sumx[oldVtx];
    for(int ii=0;ii<n;ii++){
        sumx[ii]-=m[ii][oldVtx];

        if(ii!=oldVtx&&minx[ii]==m[ii][oldVtx]){
            nmin[ii]--;
            //cout<<"nminDist="<<nminDist<<endl;
            if(x[ii]&&minx[ii]==minDist)          //The number of minimum distances is going down
                nminDist--;
            //cout<<"nmin["<<ii<<"]="<<nmin[ii]<<"x[ii]="<<x[ii]<<"oldVtx="<<oldVtx<<endl;
            if(nmin[ii]==0){                          //ii has no more minimum-value edges to our construction
                //cout<<"nminDist="<<nminDist<<endl;
                recomputeMinDistAndNMin(ii);
            }
        }
    }
    assert(nminDist>=0);
    if(nminDist==0){
        //We regenerate minDist and nMinDist
        minDist=INT_MAX;
        #ifdef DISABLE_FAST_ORDERED_LIST
        for(int ii=0;ii<n;ii++)if(x[ii])
        #else
        //std::set is using too many iterators, quite heavy
        //set<int>::iterator xIterator=selected.begin();for(int ii=*xIterator;xIterator!=selected.end();++xIterator,ii=*xIterator)
        for(int ii=selected.firstElem();ii<n;ii=selected.nextElem(ii))
        #endif
        {
            if(minx[ii]<minDist){
                minDist=minx[ii];
                nminDist=nmin[ii];
            }
            else if(minx[ii]==minDist){
                nminDist+=nmin[ii];
            }
        }
        nminDist=nminDist/2;//each vertex that belongs to an edge of minimum distance was counted twice
    }

    if(it%CHECKRATE==0){                   //from time to time, we check the streamlining calculations
        assert(checkMMDP(x,minDist));  //all asserts can however be disabled with #define NDEBUG
        assert(checkMDP(x,sumDist));
    }
}

//See Algorithm 3 for details on Algorithm below
void updateWhenAdded(int newVtx){
    int i=0;
    if(minx[newVtx]<minDist){
        minDist  = minx[newVtx];
        nminDist = nmin[newVtx];
    }
    else if(minx[newVtx]==minDist)
        nminDist += nmin[newVtx];
    for(i=0;i<n;i++)if(i!=newVtx){
            sumx[i]+=m[i][newVtx];
            if(x[i])
                sumDist+=m[i][newVtx];
            if(m[i][newVtx]<minx[i]){
                minx[i]=m[i][newVtx];
                nmin[i]=1;
            }
            else if(m[i][newVtx]==minx[i])
                nmin[i]++;
        }


}


/*
 * ==================================================================
 *
 *         ROUTINES For selecting the next vertex to add/drop
 *
 * ==================================================================
 */


int acceptableDropMove(int chosenVtx){
    return (dropTabu[chosenVtx]<it);
}
void dropVertexOldest(){
    int chosen=-1;
    for(int ii=0;ii<n;ii++)
        if(x[ii])
            if(acceptableDropMove(ii)){
//                Unless interleaved, the following should be true
//                assert(chosen==-1);
                chosen=ii;
//              The oldest one is available, it should be the only one unless interleaved
                break;
            }
    drop(chosen);
    updateWhenDropped(chosen);
}




//Returns 1 if not tabu (the tabu list mark is higher than the current iteration) or if move is best of best
int acceptableAddMove(int chosenVtx, int & asspiredStatus){
    if(addTabu[chosenVtx]<it){
        asspiredStatus = 0; //The move is really not Tabu
        return 1;
    }
    if(minDist>bestMinDistUntilNow())
        if(minx[chosenVtx]>bestMinDistUntilNow()){
            asspiredStatus = 1;
            return 1;
        }
    return 0;
}

//Attention to this criterion in the if below
#define NMIN 0
//#define NMIN (minx[i]==maxmin&&nmin[i]<nmin[chosenVtx])

int chooseAddGreedyMaxMin(){
    int i,chosenVtx=-1;
    double maxmin=INT_MIN,maxsum=INT_MIN;
    int chosenRnd=0;                                     //The last (random) minimization criterion
    int asspiredStatusValue=0;
    for(i=0;i<n;i++)
      if(x[i]==0)                                        //Only unselected vertices can be added
       if(acceptableAddMove(i,asspiredStatusValue)){     //They should be not Tabu, or at least asspired
        if(minx[i]>maxmin||NMIN){
            chosenVtx   = i;
            maxmin  = minx[i];
            maxsum  = sumx[i];
            if(asspiredStatusValue==0)
                chosenRnd = randomInt(1000);//First choosing criterion: nmin[i] (the number of points connected to i via the smallest distance), then a random number
            else
                chosenRnd = INT_MAX;       //In asspiration situation, the move is not prioritary
                                           //=> in case of equality, INT_MAX is bigger

        }
        else
            if(minx[i]==maxmin){
                if(sumx[i]>maxsum){
                    chosenVtx   = i;
                    maxsum  = sumx[i];
                    chosenRnd = randomInt(1000);;
                }
                else if(sumx[i]==maxsum){
                    int nowRand = randomInt(1000);;
                    if(nowRand<chosenRnd){
                        chosenVtx   = i;
                        //Seems useless maxsum  = sumx[i];
                        chosenRnd = nowRand;
                        //cout<<"aaa"<<endl;
                        chosenRnd = randomInt(1000);
                    }
                }
            }
      }//if(x[i]==0)
    return chosenVtx;
}

int firstPointFurthestFromAll(){
    double maxmin=INT_MIN,maxsum=INT_MIN;
    int maxrand=0;
    double currentMin;
    double currentSum;
    int currentRand;
    int bestPoint = -1; //to return this
    for(int ii=0;ii<n;ii++){
        currentMin=INT_MAX;
        currentSum=0;
        //computing currentMind and currentSum, see Constructive Algorithm
        for(int jj=0;jj<n;jj++)
            if(ii!=jj){
                currentSum+=m[ii][jj];
                if(m[ii][jj]<currentMin)
                    currentMin=m[ii][jj];
            }
        //cout<<currentMin<<endl;
        currentRand=randomInt(1000);
        if(currentMin>maxmin||
          (currentMin==maxmin&&currentSum>maxsum)||
          (currentMin==maxmin&&currentSum==maxsum&&currentRand>maxrand)
          ){
            bestPoint=ii;
            maxmin=currentMin;maxsum=currentSum;maxrand=currentRand;
        }

    }
    //cout<<"First chosen point has maxmin="<<maxmin<<" and maxsum="<<maxsum<<endl;
    return bestPoint;

}

void firstTwo(){
    int chosen;
    int chosen2 ;
    #ifdef STARTPOINT1
    chosen    = randomInt(n);
    #endif
    #ifdef STARTPOINT2
    chosen    = firstPointFurthestFromAll();
    #endif
    cout<<chosen<<endl;
    add(chosen);
    //x[chosen] = 1;
    int i;
    minx[chosen]=INT_MAX;
    sumx[chosen]=0;
    for(i=0;i<n;i++)
        if(i!=chosen)
        {
            minx[i] = m[i][chosen];
            sumx[i] = minx[i];
            nmin[i] = 1;
        }
    chosen2 = chooseAddGreedyMaxMin();
    add(chosen2);
    //x[chosen2]=1;
    assert(chosen!=chosen2);
    updateWhenAdded(chosen2);
}


void addVertexConstructive(){

    assert(h<k);
    if(h==0){
        firstTwo();
        h+=2;
    }
    int chosen=chooseAddGreedyMaxMin();
    add(chosen);
    //x[chosen]=1;
    updateWhenAdded(chosen);
}


/*
 * ==================================================================
 *
 *         ROUTINES for the drop-add LS, iterating the selections
 *
 * ==================================================================
 */

//*Simple constructive algorithm*//
//level = the number of vertices to be selected = the number of constructive steps
void constructiveAlg(int level){
    setGenericTabuTenure();     //Even if no vertex selected, we take care to have moves Tabu
    do{
        addVertexConstructive();
        h++;
    }while(h<level);
    checkNewBetterSolution();
}
void constructiveAlg(){
    constructiveAlg(k);
}

//h is the number of vertices taken into account here
//The procedure considers the h already selected points and applies local search on them
void dropAddLocalSearch(){
    genericDropTenure = h-1;//Maybe not useful (already defined in the top of the paper)
    genericAddTenure  = 0;
    //maxIter=5;
    while(itersNoNewBetter()<maxIter*n){
        dropVertexOldest();
        h=h-1;
//      cout<<h<<"vertices ("<<minDist<<")";
        addVertexConstructive();
        h=h+1;
//      cout<<h<<"vertices ("<<minDist<<")";
        checkNewBetterSolution();
    }
}
//level is the number of vertices that need to be selected
void dropAddAlg(int level){
    constructiveAlg(level);
    dropAddLocalSearch();
}
void dropAddAlg(){
    dropAddAlg(k);
}


