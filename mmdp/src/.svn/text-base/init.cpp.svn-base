#include "init.h"

void usage()
{
  std::cerr <<
    "USAGE: mmdp input_filename output_folder" << std::endl;
}
void resetIterAndTabu(){
    int i;
    for(i=0;i<n;i++){
        iter[i]=0;
        dropTabu[i]=0;
        addTabu[i]=0;
    }
    it = 1;
}
void initMatrix(){
    m = new double*[n];
    for (int ii=0;ii<n;ii++)
        m[ii]= new double[n];
    x = new int[n];
    minx = new double[n];
    sumx = new double[n];
    nmin = new int[n];
    iter = new int[n];
    dropTabu = new int[n];
    addTabu  = new int[n];
    for (int ii=0;ii<n;ii++){
        x[ii]=0;
        minx[ii]=INT_MAX;
        sumx[ii]=0;
        nmin[ii]=0;
    }
    resetIterAndTabu();
}

void printMatrix(ostream & myout, int*x){
    int i,j;
    myout<<"========================================================THE SOLUTION ON THE MATRIX============"<<endl;

    for(i=0;i<n;i++)
        if(x[i]){
            for(j=0;j<n;j++)
                if(x[i]&&x[j])
                    myout<<setw(10)<< setprecision(7)<<m[i][j];
            myout<<endl;
        }

    myout<<"================================ALL MATRIX WITH SOLUTION HIGHLIGHT========================="<<endl;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++)
            if(x[i]&&x[j])
                myout<<setw(1)<<"G"<<setw(8)<<setprecision(7)<<m[i][j]<<setw(1)<<"S";
            else
                myout<<setw(10)<< setprecision(7)<<m[i][j];
        myout<<endl;
    }
}
//*/

void readGlover(ifstream& myfile,int dimensions){
    cout<<"'Glover/Geo' matrix....";
    double**points;
    points = new double*[n];
    for (int ii=0;ii<n;ii++)
        points[ii]= new double[dimensions];
    int i,j;
    for(i=0;i<n;i++){
        int readElem;
        myfile >> readElem;
        assert(i==readElem);
        //cout<<i<<": ";
        for(int dim=0;dim<dimensions;dim++){
            myfile >> points[i][dim];
            //cout << points[i][dim]<<" ";
        }
        assert(!myfile.fail()&&!myfile.eof());
        //cout<<endl;
    }
    double squareSum;
    for(i=0;i<n;i++)
        for(j=i;j<n;j++){
            squareSum=0;
            for(int dim=0;dim<dimensions;dim++)
                squareSum+=(points[i][dim]-points[j][dim])*(points[i][dim]-points[j][dim]);
            m[i][j]=sqrt(squareSum);
            m[j][i]=m[i][j];
        }
}

void readRun(ifstream& myfile){
    cout<<" 'Run' matrix...";
    int i,j;
    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            if(i==0&&j==1)
                j++;
            //cout<<j<<endl;
            int firstInLine,seccondInLine;
            myfile>>firstInLine;
            myfile>>seccondInLine;
            assert(i==firstInLine);
            assert(j==seccondInLine);
            myfile>>m[i][j];
            assert(!myfile.fail()&&!myfile.eof());
            m[j][i]=m[i][j];
        }
        m[i][i]=0;
    }

}

//This routine is used to read a generic matrix. First, it decides which format it uses (glover or run).
void readMatrix(char* filename){
    cout<<"Reading "<<filename<<", a ";
    string line;
    ifstream myfile(filename,ios::in);
    if (myfile.is_open()){
      //Read the number of vertices, allways the first
      getline (myfile,line);
      istringstream iss(line);
      iss >> n;
      assert(iss);
      initMatrix();
      //Check (second line) if it is a file of type Glover or Run
      getline (myfile,line);
      istringstream iss2(line);
      int firstOnSecondLine;
      iss2>>firstOnSecondLine;
      //Decide whether it is glover or run
      if(firstOnSecondLine!=0)
            readGlover(myfile,firstOnSecondLine);
      else{//The line should be like: 0 1 45.0
            int secondOnSecondLine;
            int thirdOnSecondLine;
            iss2>>secondOnSecondLine;
            assert(secondOnSecondLine==1);
            iss2>>thirdOnSecondLine;
            m[0][1]=thirdOnSecondLine;
            m[1][0]=thirdOnSecondLine;
            assert(!myfile.fail()&&!myfile.eof());
            readRun(myfile);
      }
      myfile.close();
    }
    else{
      cout << "Unable to open '"<<filename<<"'"<<endl;
      exit(1);
    }
    cout<<"Done (read)"<<endl;
    //printMatrix();
}
void getInstanceNameAndK(char*filename,string&lastWord){
    string allPath(filename);
    int lastDelim=allPath.find_last_of("/");
    lastWord=allPath.substr(lastDelim+1);
    int lastSpace=lastWord.find_last_of(" ");
    int lastPoint=lastWord.find_last_of(".");
    k=atoi(lastWord.substr(lastSpace+1,lastPoint-lastSpace).c_str());
    if(k<=10)
        k=n/10;
    else
        k=3*n/10;
    lastWord=lastWord.substr(0,lastPoint);
    //lastWord is the instance, ie. it should be equal with the variable instance
    if(lastWord.find("Glover")!=string::npos)//We use the instance name(what else?)
        k=6;
}



void readConfigFile(const char* configFile){
    ifstream fConf(configFile);
    if(!fConf.is_open()) cerr<<"I cannot open the configuration file: "<<configFile<<endl;
    assert(fConf.is_open());
    string currentLine;
    while(!fConf.eof()){
        getline(fConf,currentLine);
        stringstream ss(currentLine);
        string word, wordLast, currentKey = "", currentVal = "";
        if(ss >> word)
            currentKey = currentKey + word;
        else
            continue;//Empthy line
        if(ss>>word)
            wordLast=word;
        else{
                cerr<<"Warning: ignoring "<<currentKey<<"no value for it"<<endl;
                continue;
        }
        while ( ss >> word ){
            currentKey=currentKey+" "+wordLast;
            wordLast = word;
        }
        currentVal = wordLast;
        if(configs.find(currentKey)!=configs.end())
            cerr<<"Warning, the config file has two values for (I take the last): "<<currentKey<<endl;
        configs[currentKey]=currentVal;
    }

    fConf.close();
}

template <class myType>
int getVarFromConfigTable(string myKey, myType& variable){
    map<string,string>::iterator it;
    it = configs.find(myKey);
    if(it==configs.end()){
        cout<<"Attention: NO value for "<<myKey<<" in the configuration file. The default value is "<<variable<<endl;
        return 0;
    }
    stringstream ss(it->second);
    ss>>variable;
    return 1;
}
//This function supposes that the configs map is initialized (e.g. with readConfigFile())
void readConfigParams(){
    timeMax = 0;
    getVarFromConfigTable<long>("timeMax",timeMax);
    getVarFromConfigTable<double>(instance,bestKnown);
    getVarFromConfigTable<double>("epsilon",epsilon);
    getVarFromConfigTable<int>("srandSeed",srandSeed);
    getVarFromConfigTable<int>("maxNoGain",maxNoGain);
    getVarFromConfigTable<int>("maxIter",maxIter);
}



///This routine allocates memory, reads instances and initializes other variables
void init(int argc, char**argv, string configFile){
    if(argc<3){
        usage();
        exit(EXIT_FAILURE);
    }
    readMatrix(argv[1]);

    for(int i=0;i<n;i++)
        x[i]=0;
    getInstanceNameAndK(argv[1],instance);
    outputFolder=argv[2];
    //system(("mkdir "+outputFolder+"/tex").c_str()); //LINUX
    //system(("mkdir "+outputFolder+"\\tex").c_str());//Windows
    //mkdir (("mkdir "+outputFolder+"/tex").c_str(), O_CREAT);
    texFile.open((outputFolder+"/tex/"+instance+".tex").c_str(),ios::out);
    if(texFile.fail())
		cerr<<"*********** *********** *********** ***********\nAttention: I can NOT write to the output Tex file (maybe the output folder does not exist?)"<<endl<<"*********** *********** *********** ***********\n";
    texFile<<"%Input file='"<<argv[1];
    texFile<<"(Vertices="<<n<<")\n";
    readConfigFile(configFile.c_str());
    readConfigParams();
}
