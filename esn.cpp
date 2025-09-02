#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#include <armadillo>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "dynamics.hpp"
#include "esn.hpp"
using namespace std;


typedef std::mt19937 RNG;

Esn::Esn(int size, arma::mat input, bool train){
    arma::arma_rng::set_seed_random();
    x = arma::vec(size).randu();
    arma:: mat a = arma::mat(size,1).fill(0.02);
    arma:: mat b = arma::mat(size,1).fill(0.01);
    x = x%a-b;

    arma:: sp_mat w1 = arma::sprandn(size,size,0.15); //0.6 cell dynamics //0.15 cdw

    arma::cx_vec eigval = arma::eigs_gen(w1, 1); 
    double lam = eigval(0,0).real();
    w1 = w1/abs(lam);
    arma:: mat w2 = arma::mat(w1);
    arma:: mat alpha = arma::mat(size,size).fill(0.2); //0.79 cell dynamics //0.19 cdw, 0.53 cdw new, 0.1 test //0.88
    w = w2%alpha;

    in = input;

    win = arma::mat(size,in.n_cols).randn();

    RNG rng;
    std::random_device seed;
    rng = RNG(seed());
    uniform_real_distribution<double> rd(-1,1);   
    uniform_real_distribution<double> rd2(-0.05,0.05);
    
    for(int i = 0; i<size; i++){
        arma:: vec scale = arma::vec(10).randn();
        for(int j = 0; j<49;j++){
            if (j ==0){
                win.col(j).row(i) = scale(0)*0.7; //0
            }
            else if(j== 18||j== 24 ||j== 25||j== 31){ //1
                win.col(j).row(i) = scale(1)*0.6; 
            }
            else if(j == 17 || j== 19||j== 30||j== 32){ //2
                win.col(j).row(i) = scale(2)*0.5; 
            }
            else if(j==11||j==23||j==26||j==38){ //3
                win.col(j).row(i) = scale(3)*0.3;
            }
            else if(j==10||j==12||j==16||j==20||j== 29||j== 33||j== 37||j== 39){ //4
                win.col(j).row(i) = scale(4)*0.2;
            }
            else if(j== 9||j== 13||j== 36||j== 40){ //5
                win.col(j).row(i) = scale(5)*0.16; 
            }
            else if(j== 4||j== 22||j== 27||j== 45){ //6
                win.col(j).row(i) = scale(6)*0.14; //
            }
            else if(j== 3||j== 5||j== 15||j== 21||j== 28||j== 34||j== 44||j== 46){ //7
                win.col(j).row(i) = scale(7)*0.12; 
            }
            else if(j== 2||j== 6||j== 8||j== 14||j== 35||j== 41||j== 43||j== 47){ //8
                win.col(j).row(i) = scale(8)*0.1; 
            }
            else{ //11
                win.col(j).row(i) = scale(9)*0.08;
            }
        }
    }
   
    cout<<win.row(0).n_cols<<endl;
    win.save("./esn/win.dat",arma::raw_ascii);
    w.save("./esn/w.dat",arma::raw_ascii);
}

Esn::Esn(){
    arma::arma_rng::set_seed_random();
    w.load("./esn/w.dat");
    win.load("./esn/win.dat");
    wout.load("./esn/wout.dat");
    // x.load("/home/khoan/esnDynamics/esn/x.dat");

    x = arma::mat(w.n_cols,1).randu(); 
    arma:: mat a = arma::mat(w.n_cols,1).fill(0.02);
    arma:: mat b = arma::mat(w.n_cols,1).fill(0.01);
    x = x%a-b;
    // x = arma::mat(w.n_cols,1).fill(0);
    
    
}

Esn::Esn(int L){
    arma::arma_rng::set_seed_random();
    w.load("./esn/w.dat");
    win.load("./esn/win.dat");
    wout.load("./esn/wout.dat");
    x = arma::mat(w.n_cols,1).randu(); 
    arma:: mat a = arma::mat(w.n_cols,1).fill(0.02);
    arma:: mat b = arma::mat(w.n_cols,1).fill(0.01);
    x = x%a-b;
    x_arr = new arma::mat[L*L];
    for(int i = 0; i<L*L; i++){
        x_arr[i] = arma::mat(w.n_cols,1).randu();
        x_arr[i] = x_arr[i]%a - b;
    }
    
}

void Esn::train(int steps, int numpoints){

    int len = 400;
    int size = x.n_rows;
    double a = 0.3;
    int start = 15;

    arma:: mat M = arma::mat(steps-(numpoints*start)-numpoints,size).fill(0);
    arma:: mat T = arma::mat(steps-(numpoints*start)-numpoints,1).fill(0);    
    
    RNG rng;
    std::random_device seed;
    rng = RNG(seed());
    uniform_real_distribution<double> rd(0.00001,0.001);        


    FILE* data = fopen("internal.dat","w");
    for(int k = 0; k<numpoints; k++){
        arma::arma_rng::set_seed_random();
        x = arma::vec(size).randu();
        arma:: mat a = arma::mat(size,1).fill(0.02);
        arma:: mat b = arma::mat(size,1).fill(0.01);
        x = x%a-b;
        // x = arma::mat(w.n_cols,1).fill(0);
        int init = len*k;
       
        for(int i = init; i<init+len-1; i++){
            double v = rd(rng);
            for(int j = 0; j<size; j++){
                x(j) = tanh(arma::dot(x,w.col(j))+ arma::dot(in.row(i),win.row(j)) + v);
                
                if(i >=init+start){
                    M(i-start*(k+1)-k,j) = x(j);                
                
                }
            }
            if(i>=init+start){
               T(i-start*(k+1)-k,0) = (in(i+1,0));
               
                
            }
       
            fprintf(data,"%d\t%f\n", i, x(0));
        }
    }
    cout<<"wout"<<endl;
    wout = pinv(M) * T;

    wout.save("./esn/wout.dat",arma::raw_ascii);
    x.save("./esn/x.dat",arma::raw_ascii);
}
double Esn::out(arma::vec input){   
    int size = x.n_rows;
    for(int j = 0; j<size;j++){
        x(j) = tanh(arma::dot(x,w.col(j)) + arma::dot(input,win.row(j)));
    }
    
    return (arma::dot(wout,x));
}

int* Esn::neighborESN(Index* lat, int L, int point){
    int rad = 3;
    static int nn[48];
    int xp = mod(point,L);
    int yp = point/L;
    int counter = 0;
  
    for(int y = -rad; y<=rad; y++){
        for(int x = -rad; x<=rad; x++){
            if( y==0 && x ==0){
            }
            else{
                nn[counter] = L*mod(yp + y,L) + mod(xp + x,L);
                counter++;
            }
        }
    }

    return nn;
}

int Esn::mod(int x, int m)
{
    if (x >= 0 && x < m)
        return x;

    else if (x < 0)
        return m - 1 - mod(-1 - x, m);

    else
        return x % m;
}

void Esn::latEsnUpdate(Index *lat, int L){
    int neighbors = 49;
    Index* tmp = new Index[L*L];
    memcpy(tmp, lat, L * L * sizeof(Index));
    arma::vec input = arma::vec(neighbors).fill(0);

    for(int i = 0; i<L*L;i++){
        x = x_arr[i];
        int* nn = neighborESN(lat,L,i);
        input(0) = lat[i].val;
        for(int x = 0; x<neighbors-1; x++){
            input(x+1) = lat[nn[x]].val;
        }

        tmp[i].val = out(input);
        x_arr[i] = x;
        
    }
    memcpy(lat,tmp, L * L * sizeof(Index));
    delete[] tmp;
}
