#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "nr3.h"
#include "quadrature.h"
#include "derule.h"
#include "utilities.h"


using namespace std;
double f1(double x, double delta){
    return cos(pow(x,2))*exp(-x);
}
double f2(double x, double delta){
    return sqrt(x)*cos(pow(x,2))*exp(-x);

}
double f3(double x, double delta){
    return 1000*exp(-1/x)*exp(-1/(1-x));

}
double f4(double x, double delta){
    return (1/(sqrt(x)))*cos(pow(x,2))*exp(-x);

}
int main() {
    double b= 1;
    double a = 0;
    double hmax_mild = 3.7;
    double hmax_strong = 4.3;
    // Function 1:  
    cout << "Function 1: \n";
    {
        // Instantiate DErule object.
        DErule<double (double,double)> d(f1,a,b,hmax_mild);
        int N = 10; // Iterations:
        for (int i =0; i <N;i++){
        Doub res1 =  d.next(); 
            cout << res1 << "\n";
        }
    }
    // Function 2:  
    cout << "Function 2: \n";
    {
        // Instantiate DErule object.
        DErule<double (double,double)> d(f2,a,b,hmax_mild);
        int N = 10; // Iterations:
        for (int i =0; i <N;i++){
        Doub res1 =  d.next(); 
            cout << res1 << "\n";
        }
    }
    // Function 3:  
    cout << "Function 3: \n";
    {
        // Instantiate DErule object.
        DErule<double (double,double)> d(f3,a,b,hmax_strong);
        int N = 10; // Iterations:
        for (int i =0; i <N;i++){
        Doub res1 =  d.next(); 
            cout << res1 << "\n";
        }
    }
    // Function 4:  
    cout << "Function 4: \n";
    {
        // Instantiate DErule object.
        DErule<double (double,double)> d(f4,a,b,hmax_strong);
        
        int N = 10; // Iterations:
        for (int i =0; i <N;i++){
        Doub res1 =  d.next(); 
            cout << res1 << "\n";
        }
    }
    return 0;
}
