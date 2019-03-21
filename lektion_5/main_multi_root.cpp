#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "nr3.h"
#include "utilities.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"

using namespace std;
int main() {
	VecDoub inputVector(8);	
	double pi = 3.14159265359;
	enum { L0=0,L=1,p=2,x=3,theta=4,phi=5,a=6,H=7};
	// Define problem constants:
	const double v = 120.0;
	const double k = 2.5;
	const double w = 4.0;
	const double alpha = 0.0000002;
	// Define problem paramaters
	double d = 30.0;
	double n = 500.0;
	auto func = [&] (VecDoub_IO in){
		VecDoub out(8);
		out[0]=in[a]*(cosh((in[x]/in[a])-1))-in[p];// Kædeligningen (10)
		out[1]=2.0*in[a]*sinh(in[x]/in[a])-in[L]; // Buelængden (11)
		out[2]=2.0*in[x]+2.0*k*cos(in[theta])-d;
		out[3]=in[p]+k*sin(in[theta])-n;
		out[4]=sinh(in[x]/in[a])-tan(in[phi]);
		out[5]=(1.0+v/(w*in[L0]))*tan(in[phi])-tan(in[theta]);
		out[6]=in[L0]*(1.0+alpha*in[H])-in[L];
		out[7]=(w*in[L0])/(2.0*sin(in[phi]))-in[H];
		return out;
	};
	Bool check= false;
//	auto res = func(inputVector);
	// Defines initial guesses:
	inputVector[L0] = 30;
	inputVector[L] = 30;
	inputVector[p] = 1.1;
	inputVector[x] = 15.0;
	inputVector[theta] = pi/5;
	inputVector[phi] = pi/20;
	inputVector[a] = 40;
	inputVector[H] = 5;
	try
	{
	newt(inputVector,check,func);
	}
	catch (int e)
	{
	    cout << "An exception occurred. Exception Nr. " << e << '\n';
	
	}
	//if(!check){	
	cout << "L0 ="<< inputVector[L0]<<endl;
	cout << "L ="<< inputVector[L]<<endl;
	cout << "p ="<< inputVector[p]<<endl;
	cout << "x="<< inputVector[x]<<endl;
	cout << "theta ="<< inputVector[theta]<<endl;
	cout << "phi ="<< inputVector[phi]<<endl;
	cout << "a ="<< inputVector[a]<<endl;
	cout << "H ="<< inputVector[H]<<endl;
    	//util::print(inputVector);
	//}
    	return 0;
}

