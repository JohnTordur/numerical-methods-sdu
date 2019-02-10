#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "nr3.h"
#include "ludcmp.h"
#include "svd.h"
#include "utilities.h"


using namespace std;
MatDoub makeAPolynomialFromNdegree(int n, VecDoub x){
    // First allocate the A matrice
    MatDoub *A;
    A = new MatDoub(x.size(),n); // = new MatDoub(x.size(),n);
    // Loop through all the rows:
    for(int r = 0; r < x.size();r++){
        // Generate a row at a time:
        (*A)[r][0] =1.0; 
        for( int i = 1; i < n; i++){
           (*A)[r][i] = pow(x[r],(double)i);
        }
    }
    return *A;
}
void message(string input){
    auto sbar = "----------------------------------\n";
    cout << sbar<<input << '\n'<<sbar;
}
int main() {
    // Add additional scoping to avoid unintentional errors
    // betweem the filip problem and the pontius problem.
    {
        VecDoub xPont(40); VecDoub yPont(40);
        ifstream Pont("PontiusData.dat");
        for(int i = 0; i < 40; i++) {
            Pont >> yPont[i];
            Pont >> xPont[i];
        }
        // Generate an A matrix of the form 1 x x^2 x^3 ... x^n
        //                                  1 x x^2 x^3 ... x^n
        //                                      ....    ... ...
        //                                      ....    ... ...
        //                                  1 x x^2 x^3 ... x^n
        auto A = makeAPolynomialFromNdegree(3,xPont);

        auto AT= util::Transpose(A);

        // Calculate Alpha as transpose(A)*A
        auto Alpha =AT*A;
        // Calculate Beta as transpose(A)*b
        auto Beta =AT*yPont; 

        LUdcmp ldcmp=LUdcmp(Alpha);
        VecDoub x(3);
        ldcmp.solve(Beta,x);
        util::print(x);

    }
    
    {
        VecDoub xFilip(82); VecDoub yFilip(82);
        ifstream Filip("FilipData.dat");
        for(int i = 0; i < 82; i++) {
            Filip >> yFilip[i];
            Filip >> xFilip[i];
        }
        // Generate an A matrix of the form 1 x x^2 x^3 ... x^n
        //                                  1 x x^2 x^3 ... x^n
        //                                      ....    ... ...
        //                                      ....    ... ...
        //                                  1 x x^2 x^3 ... x^n
        auto A = makeAPolynomialFromNdegree(11,xFilip);
        //util::print(A);
        auto AT= util::Transpose(A);

        // Calculate Alpha as transpose(A)*A
        auto Alpha =AT*A;
        // Calculate Beta as transpose(A)*b
        auto Beta =AT*yFilip; 
        LUdcmp ldcmp=LUdcmp(Alpha);
        VecDoub x(11);  
        ldcmp.solve(Beta,x);
       // SVD svd =SVD(Alpha);
        //svd.solve(Beta,x);
        //auto rank = svd.rank(0.00001);
        util::print(x);
    }
    return 0;
}
