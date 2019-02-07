#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "utilities.h"

using namespace std;
	
int main() {



	// Exercise 1:
	// Solve A x = b using LU decomposition, and print the result.

	MatDoub A(3,3);
	A[0][0] = 1.0;	A[0][1] = 2.0;	A[0][2] = 3.0;
	A[1][0] = 2.0;	A[1][1] = -4.0;	A[1][2] = 6.0;
	A[2][0] = 3.0;	A[2][1] = -9.0;	A[2][2] = -3.0;

	VecDoub b(3);
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	VecDoub x(3);
	util::print(b);
	LUdcmp ldcmp=LUdcmp(A);
	// evaluate x
	ldcmp.solve(b,x);
	// print x
	util::print(x);
	// Verify that A == LU
	
	return 0;
}
