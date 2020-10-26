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

    cout <<" Enter the value for sigma (>=1, Rayleigh only!): ";
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
            sum2 += (1/del)*beta(n,m,k,t*(1/del),theta,c,sig,del,base)*(1/del)*series(n,t-m);
            cout << setw(12) << nn << setw(12) << n << setw(12) << m << setw(13)<< rect(t, 0.1, 9.9)
                 << setw(18) << series(n,t) << setw(18) << sum1 << setw(18) << sum2 << setw(18) << findError(sum2,sum1);
            cout << endl;
        }
    }
    return 0;
}

/*
 *   Enter kernel: PL
 Enter the value for t: 1.2
 Enter the value for n (n>=0): 5
 Enter the value for N (N>=1): 4
 Enter the value for r: 4
 Enter the value for k: 0.8
 Enter the value for theta: 0.7
 Enter the value for c: 1.1
 Enter the value for sigma (>=1): 2
 Enter the value for delta: 0.001
           N          n         r              rect            gamma_n              h          h_n           error

           4           1    0.000000     1.000000          0.000000          0.000000          0.000000          0.000000
           4           2    0.000000     1.000000          0.800000          0.000116          0.002786          0.002670
           4           2    1.000000     1.000000          0.800000          0.000145          0.003482          0.003337
           4           3    0.000000     1.000000          0.660000          0.000241          0.005781          0.005540
           4           3    1.000000     1.000000          0.660000          0.000244          0.005850          0.005606
           4           3    2.000000     1.000000          0.660000          0.000244          0.005850          0.005606
           4           4    0.000000     1.000000          0.282667          0.000285          0.006835          0.006550
           4           4    1.000000     1.000000          0.282667          0.000285          0.006839          0.006554
           4           4    2.000000     1.000000          0.282667          0.000285          0.006839          0.006554
           4           4    3.000000     1.000000          0.282667          0.000285          0.006839          0.006554

Process finished with exit code 0

 */