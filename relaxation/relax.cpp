#include "nr3.h"
#include "tridag.h"

#include "utilities.h"

Doub F(Doub ydot, Doub y, Doub x) { return -cos(y) * sin(ydot); }
Doub Fy(Doub ydot, Doub y, Doub x) {
    // Derivative with respect to y
    // sin(y')*sin(y)
    return sin(ydot) * sin(y);
}
Doub Fyy(Doub ydot, Doub y, Doub x) {
    return -cos(y) * cos(ydot);
    // Derivative with respect to y'
    //-cos(y)*cos(y')
}

// void calculateJacobian(VecDoub_I y,VecDoub_I avec, VecDoub_I bvec, VecDoub_I
// cvec,

int main() {
    Int N = 100;
    Doub a = 0.0;
    Doub b = 10.0;
    Doub h = (b - a) / (Doub)N;
    // Calculate Alpha
    Doub alpha = 0;
    // Calculate Beta.
    Doub beta = 3;
    // Compute straight line approximation(initial guess):
    VecDoub yguess(N, 0.0);
    // Compute alle x's
    VecDoub x(N, 0.0);
    for (int i = 0; i < N; i++) {
        x[i] = a + i * h;
    }
    // util::print(x);
    yguess[0] = alpha;
    for (int i = 1; i < N - 1; i++) {
        yguess[i] = alpha + (i / ((Doub)N)) * (beta - alpha);
    }
    yguess[N - 1] = beta;
    // util::print(yguess);
    // Set up jacobian elements:

    // vector a = i,i-1
    VecDoub avec(N, 0.0);
    avec[N - 2] = -666;
    // vector b = i,i
    VecDoub bvec(N, 0.0);
    // vector c = i,i+1
    VecDoub cvec(N, 0.0);
    // Husk at a[0] er udefineret
    for (int j = 0; j < 10; j++) {
        Doub y2 = yguess[2];
        Doub y1 = yguess[1];

        // J1,1
        bvec[0] = 2 + h * h * Fy((y2 - alpha) / 2 * h, y1, x[1]);
        // J1,2
        cvec[0] = -1 + (h / 2) * Fyy((y2 - alpha) / 2 * h, y1, x[1]);

        // JN-1,N-2
        avec[N - 1] = -1 - (h / 2) * Fyy((beta - yguess[N - 2]) / (2 * h),
                                         yguess[N - 1], x[N - 1]);
        // JN-1,N-1
        bvec[N - 1] = 2 + h * h *
                              Fy((beta - yguess[N - 2]) / (2 * h),
                                 yguess[N - 1], x[N - 1]);
        // cout << "bvec should be defined at 0 and N-1\n";
        // util::print(bvec);
        // cout << "\n";
        // cout << "avec should be defined at N-1 only\n";
        // util::print(avec);
        // cout << "\n";
        // cout << "cvec should be defined at 0 only\n";
        // util::print(cvec);
        // cout << "\n";
        for (int i = 1; i < N - 1; i++) {
            Doub y1 = yguess[i + 1];
            Doub y_1 = yguess[i - 1];
            Doub yi = yguess[i];
            Doub xi = x[i];
            // Ji,i-1
            avec[i] = -1.0 - (h / 2.0) * Fyy((y1 - y_1) / (2 * h), yi, xi);
            // Ji,i
            bvec[i] = 2.0 + h * h * Fy((y1 - y_1) / (2 * h), yi, xi);
            // Ji,i+1
            cvec[i] = -1.0 + (h / 2.0) * Fyy((y1 - y_1) / (2 * h), yi, xi);
        }

        // cout << "bvec should be defined everywhere\n";
        // util::print(bvec);
        // cout << "\n";
        // cout << "avec should be defined at 1 .. to N-1 only\n";
        // util::print(avec);
        // cout << "\n";
        // cout << "cvec should be defined at 0 to N-2\n";
        // util::print(cvec);
        // cout << "\n";

        // Define phi Vector
        VecDoub phiVec(N, 0.0);
        phiVec[1] =
            -2 * y1 + y2 - h * h * F((y2 - alpha) / (2 * h), y1, x[1]) + alpha;
        phiVec[N - 1] =
            yguess[N - 2] - 2 * yguess[N - 1] -
            h * h *
                F((beta - yguess[N - 2]) / (2 * h), yguess[N - 1], x[N - 1]) +
            beta;
        for (int i = 2; i < N - 1; i++) {
            phiVec[i] = yguess[i - 1] - 2 * yguess[i] + yguess[i + 1] -
                        h * h *
                            F((yguess[i + 1] - yguess[i - 1]) / (2 * h),
                              yguess[i], x[i]);
        }
        // Compute
        //
        for (int i = 0; i < N; i++) {
            cout << x[i] << " ";
            cout << yguess[i];
            cout << "\n";
        }
        VecDoub uvec(N, 0.0);
        // Solve
        tridag(avec, bvec, cvec, phiVec, uvec);
        // Update
        //
        // Plot:
        for (int i = 0; i < N; i++) {
            yguess[i] = yguess[i] + uvec[i];
        }

        // util::print(yguess);
    }
    return 0;
}

