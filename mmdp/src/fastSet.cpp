class fastSet{
    private:
    int * next;
    int * prev;
    int first;
    int last;
    int size;
    int empthySet ;       //If this is 1, then the set has no elements

    void initSet(){
        size = n;
        if(next==NULL){
            next = new int[size];
            prev = new int[size];
        }
        first = size;
        last = -1;
        for(int ii=0;ii<size;ii++)
            next[ii]=-1;      //next[ii]=-1 means unselectected
        empthySet = 0;
    }

    public:
    fastSet(): next (NULL), prev (NULL), empthySet(1) {};
    void reInitialize(){
//        delete next;
//        delete prev;
//        next = NULL;
        empthySet = 1;
    }
    //I suppose there is a left, not -1;
    int findLeft(int current){
        int left=current-1;
        for(;left>=0;left--)
            if(next[left]>=0)
                return left;
        cerr<<"Could not find the left one\n";
        exit(1);
    }
    //Attention: never not delete the last and only element
    void insert(int elem){
        if(empthySet){
            initSet();
            first = elem;
            last  = elem;
            prev[elem]=-1;
            next[elem]=size;//Pointing out of range
            return;
        }
        if(elem<first){     //Insertion in front of the ordered list
            prev[first]= elem;
            next[elem] = first;
            prev[elem] = -1;
            first      = elem;
        }
        else if(last<elem){  //Insertion at the end of the ordered list
            prev[elem] = last;
            next[elem] = size;
            next[last] = elem;
            last       = elem;
        }
        else{                //Insertion between two elements left and right=next[left]
            int left;
            left = findLeft(elem);
            prev [next[left]] = elem;
            next [elem]  = next[left];
            prev [elem]  = left;
            next [left]  = elem;
        }
    }
    void erase(int elemOut){
        if(elemOut == first){
            first = next[first];
            prev[first] = -1;
            next[elemOut] = -1;


        }
        else if(elemOut==last){
            next[last]        = -1;
            next[prev[last]]  = size;
            last = prev[last];
        }
        else{
            prev[next[elemOut]] = prev[elemOut];
            next[prev[elemOut]] = next[elemOut];
            next[elemOut]       = -1;      //No longer selected
        }
    }
    int nextElem(int current){
        return next[current];
    }
    int firstElem(){
        return first;
    }
};
