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
            sum2 += (1/del)*beta(n,m,k,t*(1/del),theta,c,sig,del,base)* series(n,t-m);
            cout << setw(12) << nn << setw(12) << n << setw(12) << m << setw(13)<< rect(t, 0.1, 9.9)
                 << setw(18) << series(n,t) << setw(18) << sum1 << setw(18) << sum2 << setw(18) << findError(sum2,sum1);
            cout << endl;
        }
    }
    return 0;
}

/*
 *  Enter kernel: PL
 Enter the value for t: 1.5
 Enter the value for n (n>=0): 6
 Enter the value for N (N>=1): 5
 Enter the value for r: 5
 Enter the value for k: 0.8
 Enter the value for theta: 0.7
 Enter the value for c: 1.1
 Enter the value for sigma (>=1): 2.2
 Enter the value for delta: 0.001
           N          n         r              rect            gamma_n              h          h_n           error

           5           1    0.000000     1.000000          0.000000          0.000000          0.000000          0.000000
           5           2    0.000000     1.000000          0.500000          0.000059          0.000001          0.000058
           5           2    1.000000     1.000000          0.500000          0.000118          0.000002          0.000116
           5           3    0.000000     1.000000          0.750000          0.000206          0.000004          0.000202
           5           3    1.000000     1.000000          0.750000          0.000221          0.000004          0.000217
           5           3    2.000000     1.000000          0.750000          0.000221          0.000004          0.000217
           5           4    0.000000     1.000000          0.479167          0.000278          0.000006          0.000272
           5           4    1.000000     1.000000          0.479167          0.000280          0.000006          0.000274
           5           4    2.000000     1.000000          0.479167          0.000280          0.000006          0.000274
           5           4    3.000000     1.000000          0.479167          0.000280          0.000006          0.000274
           5           5    0.000000     1.000000          0.197917          0.000303          0.000006          0.000297
           5           5    1.000000     1.000000          0.197917          0.000304          0.000006          0.000298
           5           5    2.000000     1.000000          0.197917          0.000304          0.000006          0.000298
           5           5    3.000000     1.000000          0.197917          0.000304          0.000006          0.000298
           5           5    4.000000     1.000000          0.197917          0.000304          0.000006          0.000298

 */