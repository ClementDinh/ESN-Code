#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#include <armadillo>
#include "dynamics.hpp"
#include "esn.hpp"
#include <iostream>
#include <string>
#include<bits/stdc++.h>
using namespace std;

typedef std::mt19937 RNG;

int main(int argc, char *argv[]){
    bool train = false;
    string arg3(argv[1]);
    if(arg3 == "train"){
        train = true;
    }

    int size = 40;
    int len = 400;
    double points[len][size*size];

    int numpoints = 150;
    RNG rng;
    std::random_device seed;
    rng = RNG(seed());
    int randindex[numpoints];
    int validpoints[size*size];
    
    uniform_real_distribution<double> rd2(0,size*size);

    int start = 0;

    int neighbors = 49;
    
    arma::mat input = arma::mat(len*numpoints,neighbors);
    

    Index* lat = new Index[size*size];
    fillIndex(lat,size);
    Index* esnLat = new Index[size*size];
    fillIndex(esnLat,size);

    int run = strtol(argv[3],NULL,0);
    
    string folder = "./CDW_paper_review_3_26/g_0.9/";
    string filePath = "./CDW_paper_review_3_26/g_0.9/CDW" + to_string(run) + ".dat";
 

    string myText;
    ifstream MyReadFile(filePath);
    string myText2;
    ifstream MyReadFile2(filePath);

        
    string word;
    string word2;
    for(int t = 0; t<len;t++){
        getline(MyReadFile2,myText2);
        int n = myText2.length();
        int counter = 0;
        for(int i = 0; i<n;i++){
            if(myText2[i] == ' ' || i==(n-1)){
                points[t][counter] = stod(word2);
                counter++;
                word2 = "";
            }
            else{
                word2 += myText2[i];
            }
        }
    }
    int nums =0;
    for(int n = 0; n<size*size; n++){
        double max = -10000;
        double min = 10000;
        for(int t = 0; t<len; t++){
            if(points[t][n] > max){
                max = points[t][n];
            }
            if(points[t][n] < min){
                min = points[t][n];
            }
        }

        if(abs(max-min) > 1.2){ //0.6 145 pts      
            validpoints[nums] = n;
            nums++;
        }
        
    }

    uniform_real_distribution<double> rd1(0,nums); //1.5 172, 1.4 367 176
    cout<<"validponts "<<nums<<endl;
  
        
    for(int i = 0; i<numpoints;i++){
        int num = (int)rd1(rng);
        int point = validpoints[num];
        while(point == -1){
            num = (int)rd1(rng);
            point = validpoints[num];
        }
        randindex[i] = point;
        validpoints[num] = -1;
    }
    
    for(int t = 0; t<400;t++){
        getline(MyReadFile,myText);
        int n = myText.length();
        int counter = 0;
        for(int i = 0; i<n;i++){
            if(myText[i] == ' ' || i==(n-1)){
                lat[counter].val = stod(word);
                counter++;
                word = "";
            }
            else{
                word += myText[i];
            }
        }
        
        if(t >= start && t <len+start){
            for(int k = 0; k<numpoints;k++){
                int point = randindex[k];
                int* nn = neighbor(lat,size,point);
                int ind = t+(len*k-start);
                input(ind,0) = lat[point].val; 
                for(int x = 0; x<neighbors-1;x++){
                    input(ind,x+1) = lat[nn[x]].val;
                }
            }
            
        }
        if(t==15){
            memcpy(esnLat,lat,size*size*sizeof(Index));
        }   
      
        string str = folder + "actual/CDW_time" + to_string(t) + ".dat";
        char name[100];
        strcpy(name, str.c_str());
        FILE* data = fopen(name,"w");
        writeLat(lat,size,data);
    }
    

    Esn my_esn(size);
    
    if( train == true){
        my_esn = Esn(200, input, true);
        my_esn.train(len*numpoints, numpoints);
    }

    input.save("./esn/input.dat",arma::raw_ascii);
    arma::vec sample = input.col(0);
    sample.save("./esn/sample.dat",arma::raw_ascii);
  
    bool test = false;
    string arg4(argv[2]);
    if(arg4 == "test"){
        test = true;
    }
    
    FILE* data5 = fopen("neighborstest.dat","w");
    
    if(test == true){
        time_t time1 = time(nullptr);
        for (int n=15;n<500+1;n++){
            string str = folder + "predicted/CDW_time" + to_string(n+1) + ".dat";
            char name[100];
            strcpy(name, str.c_str());
            double test = 0;
            
            if(n>0){
                fprintf(data5, "%d\t%f\t%f\n",n,esnLat[randindex[0]].val,sample[n-0]);
            }  

            FILE* data2 = fopen(name, "w");
            my_esn.latEsnUpdate(esnLat,size);       
               
            writeLat(esnLat,size,data2);
        }
        cout<<time(nullptr)-time1<<endl;
    }
 
    return 0;
}



