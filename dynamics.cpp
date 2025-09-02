#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdio>
#include "dynamics.hpp"

using namespace std;



int mod(int x, int m)
{
    if (x >= 0 && x < m)
        return x;

    else if (x < 0)
        return m - 1 - mod(-1 - x, m);

    else
        return x % m;
}

inline int index(int x, int y, int L) {
    return y * L + x;
}

void fillIndex(Index *lat, int L)
{
    for(int i = 0; i<L*L; i++){
        lat[i].x = mod(i,L);
        lat[i].y = i/L;
        lat[i].left = lat[i].y * L + mod(lat[i].x - 1, L);
        lat[i].right = lat[i].y * L + mod(lat[i].x + 1, L); 
        lat[i].bottom = mod(lat[i].y - 1, L) * L + lat[i].x;
        lat[i].top = mod(lat[i].y + 1, L) * L + lat[i].x;
        lat[i].topLeft = mod(lat[i].y + 1, L) * L + mod(lat[i].x - 1, L);
        lat[i].topRight = mod(lat[i].y + 1, L) * L + mod(lat[i].x + 1, L);
        lat[i].bottomLeft = mod(lat[i].y - 1, L) * L + mod(lat[i].x - 1, L);
        lat[i].bottomRight = mod(lat[i].y - 1, L) * L + mod(lat[i].x + 1, L);
    }

}

void writeLat(Index *lat, int L, FILE *data){
    for(int i = 0; i<L*L; i++){
        fprintf(data,"%d\t%d\t%f\t\n", lat[i].x, lat[i].y, lat[i].val);
        if(lat[i].x == L-1){
            fprintf(data,"\n");
        }
    }
    
}

void updateLat(Index *lat, int L){
    double A = 1.3;
    double D = 0.5;

    Index* tmp = new Index[L*L];
    memcpy(tmp, lat, L * L * sizeof(Index));

    for(int i = 0; i<L*L; i++){
        int left = lat[i].left;
        int right = lat[i].right;
        int top = lat[i].top;
        int bottom = lat[i].bottom;
        int topLeft = lat[i].topLeft;
        int topRight = lat[i].topRight;
        int bottomLeft = lat[i].bottomLeft;
        int bottomRight = lat[i].bottomRight;

        double sum = (lat[left].val + lat[right].val + lat[top].val + lat[bottom].val)/6.0 + 
                     (lat[topLeft].val + lat[topRight].val + lat[bottomLeft].val + lat[bottomRight].val)/12.0;
        
        tmp[i].val = A * tanh(lat[i].val) + D * (sum - lat[i].val);
    }
    memcpy(lat,tmp, L * L * sizeof(Index));
    delete[] tmp;
}

double corr(int x, Index* lat,int L){
    double avg1 = 0;
    double avg2 = 0; 
    double avg3 = 0;

    for (int i = 0; i < L * L; i++)
    {
        int x1 = i%L;
        int y = i/L;
        avg1 += lat[i].val;
        avg2 += 0.5*(lat[i].val * lat[(x1+x)%L+y*L].val + lat[i].val * lat[x1+((y+x)%L)*L].val);
        avg3 += lat[i].val*lat[i].val;

    }
    avg1 = pow(avg1 / (double)(L * L), 2);
    avg2 = avg2 / (double)(L * L);
    avg3 = avg3 / (double)(L*L);

    // if(x==0){
    //     cout<<(avg2-avg1)/(avg3-avg1)<<" ";
    //     cout<<(avg3)<<" "<<avg2<<endl;
        
    // }
    // // (avg2-avg1)/(avg3-avg1)
    // return (avg2-avg1)/(avg3-avg1);
    return (avg2);
}

double corr2(int x, Index* lat,int L){ //avg over lattice
    double avg1 = 0;
 

    for (int i = 0; i < L * L; i++)
    {
        avg1 += lat[i].val;

    }
    avg1 = avg1 / (double)(L * L);

    return abs(avg1);
}

int* neighbor(Index* lat, int L, int point){
    int rad = 3;
    static int nn[48];
    int xp = mod(point,L);
    int yp = point/L;
    // cout<<mod(yp+-2,L)<<endl;
    int counter = 0;
    // nn[0] = L*mod(yp + 1,L) + xp;
    // nn[1] = L*mod(yp -1,L) + xp;
    // nn[2] = L*yp + mod(xp + 1,L);
    // nn[3] = L*yp + mod(xp -1,L);
    


    for(int y = -rad; y<=rad; y++){
        for(int x = -rad; x<=rad; x++){
            if(y==0 && x== 0){
            }
            else{
                nn[counter] = L*mod(yp + y,L) + mod(xp + x,L);
                counter++;
            }
            // nn[counter] = L*mod(yp + y,L) + mod(xp + x,L);
            // counter++;

        }
    }
    
    return nn;
}


