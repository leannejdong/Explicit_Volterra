#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

int main() {
    using volterra::rect;
    using volterra::beta;
    using volterra::series;
    using volterra::findError;
    using volterra::norms;

    int n, nn, r;
    double k, t, theta, c, sig, del;
    std::string base;


    cout <<" Enter kernel: ";
    getline (cin, base);
    //cin >> base;
    cout <<" Enter the value for t: ";
    cin >> t;
    cout <<" Enter the value for n (n>=0): ";
    cin >> n;
    cout <<" Enter the value for N (N>=1): ";
    cin >> nn;
    cout << " Enter the value for r: ";
    cin >> r;
    cout <<" Enter the value for k: ";
    cin >> k;
    cout <<" Enter the value for theta: ";
    cin >> theta;
    cout <<" Enter the value for c: ";
    cin >> c;

    cout <<" Enter the value for sigma (>=1): ";
    cin >> sig;

    cout <<" Enter the value for delta: ";
    cin >> del;

//    double* x = new double[n];
//    x[0]=-5.1; x[1]=2.3; x[2]=3.7; x[3]=1.1; x[4]=0.7;
//    double norm1;
//    double norminf;
//    norms(x,n,norm1,norminf);
//    cout << "norm1 = " << norm1 << endl;
//    cout << "norminf = " << norminf << endl;
//    delete[] x;


    cout  << "           N          n         r              rect            gamma_n              h          h_n           error";
    cout << "  " << "  " << endl << endl;

    cout << fixed << setprecision(6);

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int n=1; n<nn+1; n++){
        for (double m=0; m<n; m++){
            sum1 += beta(n,m,k,t,theta,c,sig,del,base)* series(n,t-m);
            sum2 += (1/del)*beta(n,m,k,t*(1/del),theta,c,sig,del,base)* series(n,t/del-m);
            cout << setw(12) << nn << setw(12) << n << setw(12) << m << setw(13)<< rect(t, 0.1, 9.9)
                 << setw(18) << series(n,t) << setw(18) << sum1 << setw(18) << sum2 << setw(18) << findError(sum2,sum1);
            cout << endl;
        }
    }
    return 0;
}

/*
 *   Enter kernel: PL
 Enter the value for t: 0.4
 Enter the value for n (n>=0): 5
 Enter the value for N (N>=1): 6
 Enter the value for r: 4
 Enter the value for k: 0.6
 Enter the value for theta: 0.7
 Enter the value for c: 1
 Enter the value for sigma (>=1): 2.2
 Enter the value for delta: 0.2
           N          n         r              rect            gamma_n              h          h_n           error

           6           1    0.000000     1.000000          1.000000          0.047409          0.000000          0.047409
           6           2    0.000000     1.000000          0.400000          0.066373          0.000000          0.066373
           6           2    1.000000     1.000000          0.400000          0.066373          0.064885          0.001488
           6           3    0.000000     1.000000          0.080000          0.070165          0.097327          0.027162
           6           3    1.000000     1.000000          0.080000          0.070165          0.129770          0.059604
           6           3    2.000000     1.000000          0.080000          0.070165          0.129770          0.059604
           6           4    0.000000     1.000000          0.010667          0.070671          0.173026          0.102355
           6           4    1.000000     1.000000          0.010667          0.070671          0.183840          0.113169
           6           4    2.000000     1.000000          0.010667          0.070671          0.183840          0.113169
           6           4    3.000000     1.000000          0.010667          0.070671          0.183840          0.113169
           6           5    0.000000     1.000000          0.001067          0.070722          0.213579          0.142858
           6           5    1.000000     1.000000          0.001067          0.070722          0.216283          0.145561
           6           5    2.000000     1.000000          0.001067          0.070722          0.216283          0.145561
           6           5    3.000000     1.000000          0.001067          0.070722          0.216283          0.145561
           6           5    4.000000     1.000000          0.001067          0.070722          0.216283          0.145561
           6           6    0.000000     1.000000          0.000085          0.070726          0.230341          0.159615
           6           6    1.000000     1.000000          0.000085          0.070726          0.230882          0.160156
           6           6    2.000000     1.000000          0.000085          0.070726          0.230882          0.160156
           6           6    3.000000     1.000000          0.000085          0.070726          0.230882          0.160156
           6           6    4.000000     1.000000          0.000085          0.070726          0.230882          0.160156
           6           6    5.000000     1.000000          0.000085          0.070726          0.230882          0.160156

Process finished with exit code 0


 */