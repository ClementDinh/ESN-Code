#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#pragma once


#include <armadillo>
#include "dynamics.hpp"
class Esn{
    public:
        arma::vec x;
        arma::mat w;
        arma::mat in;
        arma::mat win;
        arma::mat wout;
        arma::mat* x_arr;

        Esn(int size, arma::mat input, bool train);
        Esn();
        Esn(int L);
        void train(int steps, int numpoints);
        double out(arma::vec input);
        void latEsnUpdate(Index *lat, int L);
        int* neighborESN(Index *lat, int L, int point);
        int mod(int x, int m);
};
