#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "nr3.h"
#include "utilities.h"


using namespace std;
MatDoub makeAPolynomialFromNdegree(int n, VecDoub x){
    // First allocate the A matrice
    MatDoub *A;
    A = new MatDoub(x.size(),n); // = new MatDoub(x.size(),n);
    // Loop through all the rows:
    for(int r = 0; r < x.size();r++){
        // Generate a row at a time:
        for( int i = 0; i < n; i++){
           (*A)[r][i] = pow(x[r],i);
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
        auto A = makeAPolynomialFromNdegree(3,xFilip);

        message("A");
        util::print(A);
        auto ATranspose = util::Transpose(A);
        message("A Tranposed:");
        util::print(ATranspose);
        cout << ATranspose.nrows()<< "\n";

        // Calculate Alpha as transpose(A)*A
        auto Alpha =A*ATranspose;
        message("Alpha:");
        util::print(Alpha);
        // Calculate Beta as transpose(A)*b
        auto Beta = ATranspose*yFilip; 
        message("Beta:");
        util::print(Beta);
    }
    VecDoub xPont(40); VecDoub yPont(40);
    ifstream Pont("PontiusData.dat");
    for(int i = 0; i < 40; i++) {
        Pont >> yPont[i];
        Pont >> xPont[i];
    }


    return 0;
}
