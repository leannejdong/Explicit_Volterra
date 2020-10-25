#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

int main() {
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

    double* x = new double[n];
    x[0]=-5.1; x[1]=2.3; x[2]=3.7; x[3]=1.1; x[4]=0.7;
    double norm1;
    double norminf;
    norms(x,n,norm1,norminf);
    cout << "norm1 = " << norm1 << endl;
    cout << "norminf = " << norminf << endl;
    delete[] x;


    cout  << "           N          n       beta_n     gamma_n      h_n          h           error";
    cout << "  " << "  " << endl << endl;

    cout << fixed << setprecision(6);

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int n=1; n<nn+1; n++){
        for (int r=0; r<n; r++){
            sum1 += beta(n,r,k,t,theta,c,sig,del,base)* series(n,t-r);
            sum2 += (1/del)*beta(n,r,k,t,theta,c,sig,del,base)* series(n,t*(1/del)-r);
            cout << setw(12) << n << setw(12) << r << setw(13) << beta(n,r,k,t,theta,c,sig,del,base) << setw(12) << series(n,t-r) << setw(13)<< sum1 << setw(13) << sum2 << setw(13) << findError(sum2,sum1);
            cout << endl;
        }
    }
    return 0;
}
