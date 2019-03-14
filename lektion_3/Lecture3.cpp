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
        for( int i = 0; i < n; i++){
           (*A)[r][i] = pow(x[r],(double)i);
        }
    }
    return *A;

}
void vecSub(VecDoub &l,VecDoub &r, VecDoub & result){

        // Euclidean norm.
        double accumSum = 0;
        for( int i = 0; i < l.size();i++){
            result[i] = l[i] - r[i];
        }
}
double vecLen(VecDoub &x){

        // Euclidean norm.
        double accumSum = 0;
        for( int i = 0; i < x.size();i++){
            accumSum += x[i] * x[i];
        }
        double norm = sqrt(accumSum);
        return norm;
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

        SVD svddcmp(A);
        VecDoub x(3) ;
        svddcmp.solve(yPont,x);
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

        SVD svddcmp(A);
        VecDoub x(11) ;
        svddcmp.solve(yFilip,x,svddcmp.eps);
        util::print(x);
        auto res = A*x;
        std::cout <<"------------------";
        util::print(res);
        std::cout <<"------------------";
        util::print(yFilip);

        std::cout <<"Residual error:";
        // Euclidean norm.
        VecDoub residual(x.size());
        vecSub(res,yFilip,residual);
        double up = vecLen(residual);   
        double low = vecLen(yFilip);
        double relativeError = up/low;
        std::cout << relativeError << endl;
        
        std::cout <<"Random model error:";
        double n = double(A.ncols());

        double  m = double(A.nrows());
        double randomModel = sqrt((m-n)/m);
        std::cout << randomModel;
    }
    return 0;
}
