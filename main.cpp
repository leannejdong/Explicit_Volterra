#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cassert>

using namespace std;

static double
kernel(
        //int,
        double k,
        double t,
        double theta,
        double c,
        double sig,
       // double,
        std:: string& base
)
{
    if(base== "PL"){
        // return k*theta*pow(c,theta)/(pow((c+(r+1)*del),(1+theta)));
        return k*theta*pow(c,theta)/(pow((c+t),(1+theta)));
    }
    else if (base=="Exp"){
        return k*theta*exp(-theta*t);
    }
    else if (base=="RL"){
        return (k*t/(pow(sig,2)))*exp(-t*t/(2*pow(sig,2)));
    }

    assert(false); // shouldn't happen
}


int main() {
    using volterra::rect;
    using volterra::beta;
    using volterra::gamma_n;
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

  //  cout <<" Enter the value for sigma (Enter a number >=1, for Rayleigh only! Otherwise, enter any key to continue): ";
  //  cin >> sig;

    cout <<" Enter the value for delta: ";
    cin >> del;

    double x[] = {1,-2,3};
    double norm1;
    double norminf;
    norms(x, 3, norm1, norminf);
    assert(x[0] == 1);
    assert(x[1] == 2);
    assert(x[2] == 3);
    assert(norm1 == 1+2+3);
    assert(norminf == 3);

    cout  << "           N          n         r              rect            gamma_n              h          h_n           error";
    cout << "  " << "  " << endl << endl;

    cout << fixed << setprecision(6);

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int n=1; n<nn+1; ++n){
        for (double m=0; m<n-1; ++m){
            auto a =
                    [&](int r, double t){
                        return kernel(k, t, theta, c, sig, base);
                    };
            sum1 += beta(n, m, a, t)* gamma_n(t-m,n);
            //sum2 += (1.0/del)*beta(n,m,a,t)*gamma_n(t*(1.0/del)-m, n);
            for (double mk = 0; mk < n; ++mk){
                sum2 += beta(n, m, a, t)*gamma_n((1.0/del)*t-m-mk, n);
            }
         //   }
            cout << setw(12) << nn << setw(12) << n << setw(12) << m << setw(13)<< rect(t, 0.1, 9.9)
                 << setw(18) << gamma_n(t, n) << setw(18) << sum1 << setw(18) << sum2 << setw(18) << findError(sum2,sum1);
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