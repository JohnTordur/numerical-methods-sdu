#include <cmath>
//#include <cstdio>
//#include <fstream>
//#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "stepper.h"
#include "odeint.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
#include "stepperdopr853.h"
#include "shoot.h"
#include "utilities.h"

using namespace std;
struct Rhs {
    Rhs() {}
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0] = y[1];
        dydx[1] = -cos(y[0]) * sin(y[1]);
    }
};

struct Load {
    Load() {}
    VecDoub operator() (const Doub x1, VecDoub_I &v) {
        VecDoub ystart(2);
        ystart[0] = 0;
        ystart[1] = v[0];
        return ystart;
    }
};

struct Score {
    Doub yt;
    Score(Doub yb) : yt(yb){}
    VecDoub operator() (const Doub xf, VecDoub_I &y) 
    {
        VecDoub error(1);
        error[0] = y[0] - yt;
        return error;
    }
};
int main() {
    Int Neq = 2;
    Doub x1 = 0.0;
    Doub x2 = 10.0;
    int N2 = 2;
    
    Load ld = Load();
    Score s= Score(3.0);
    Rhs rhs = Rhs();
    //Score s(3.0);
    VecDoub v(N2);

    v[0] = 5;
    Shoot<Load, Rhs, Score> sht(Neq, x1, x2, ld, rhs, s);
    bool check;
    newt(v,check,sht);
    if(!check){
        cout << v[0]<<"\n";
    }else{
        cout << "Shoot failed"<<endl;
    }
    // Verify:
    
    Doub h1=(x2-x1)/100.0;
    VecDoub y = ld(0.0,v);
    Doub atol(1.0e-14);
    Doub rtol(atol);
    Output out(100);
    Doub hmin = 0.0;
    Odeint<StepperDopr853<Rhs> > integ(y,x1,x2,atol,rtol,h1,hmin,out,rhs);
    integ.integrate();
   // cout << y[0]<< " Should be equal to 3 \n";
    VecDoub yres(out.count,out.ysave[0]);
   // cout << out.count << "\n";
    for(int i =0;i < out.count;i++){
    cout << out.xsave[i] << " ";
    cout << yres[i];
    cout <<"\n";
    }
    //util::print(yres);
    
    return 0;
}
