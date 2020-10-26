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


    cout  << "           N          n         r          beta_n          rect            gamma_n              h          h_n           error";
    cout << "  " << "  " << endl << endl;

    cout << fixed << setprecision(6);

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int n=1; n<nn+1; n++){
        for (double m=0; m<n; m++){
            sum1 += beta(n,m,k,t,theta,c,sig,del,base)* series(n,t-m);
            sum2 += (1/del)*beta(n,m,k,t*(1/del),theta,c,sig,del,base)* series(n,t*(1/del)-m);
            cout << setw(12) << nn << setw(12) << n << setw(12) << m << setw(13) << beta(n,m,k,t,theta,c,sig,del,base) << setw(13)<< rect(t, 0.1, 9.9)
                 << setw(18) << series(n,t) << setw(18) << sum1 << setw(18) << sum2 << setw(18) << findError(sum2,sum1);
            cout << endl;
        }
    }
    return 0;
}

/*
 *  Enter kernel: Exp
 Enter the value for t: 0.9
 Enter the value for n (n>=0): 4
 Enter the value for N (N>=1): 5
 Enter the value for r: 3
 Enter the value for k: 0.7
 Enter the value for theta: 0.8
 Enter the value for c: 1
 Enter the value for sigma (>=1): 1.3
 Enter the value for delta: 0.01
           N          n         r          beta_n          rect            gamma_n              h          h_n           error

           5           1    0.000000     0.002726     1.000000          1.000000          0.002726          0.000000          0.002726
           5           2    0.000000     0.002726     1.000000          0.900000          0.005179          0.000000          0.005179
           5           2    1.000000     0.002726     1.000000          0.900000          0.005179          0.000000          0.005179
           5           3    0.000000     0.002726     1.000000          0.405000          0.006283          0.000000          0.006283
           5           3    1.000000     0.002726     1.000000          0.405000          0.006283          0.000000          0.006283
           5           3    2.000000     0.002726     1.000000          0.405000          0.006283          0.000000          0.006283
           5           4    0.000000     0.002726     1.000000          0.121500          0.006614          0.000000          0.006614
           5           4    1.000000     0.002726     1.000000          0.121500          0.006614          0.000000          0.006614
           5           4    2.000000     0.002726     1.000000          0.121500          0.006614          0.000000          0.006614
           5           4    3.000000     0.002726     1.000000          0.121500          0.006614          0.000000          0.006614
           5           5    0.000000     0.002726     1.000000          0.027338          0.006689          0.000000          0.006689
           5           5    1.000000     0.002726     1.000000          0.027338          0.006689         -0.000000          0.006689
           5           5    2.000000     0.002726     1.000000          0.027338          0.006689         -0.000000          0.006689
           5           5    3.000000     0.002726     1.000000          0.027338          0.006689          0.000000          0.006689
           5           5    4.000000     0.002726     1.000000          0.027338          0.006689         -0.000000          0.006689

Process finished with exit code 0

 */