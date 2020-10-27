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

    int n, nn;
    double k, theta, c, sig;
    std::string base;

    base = "PL"; n = 80; nn = 5; k = 80; theta = 0.8; c = 1.1;

    cout << fixed << setprecision(6);

    int divs = 100;
    double low = -1;
    double high = n+1;
    cerr << "gamma_n(1200)=" << gamma_n(1200,n) << "\n";
    for (int i = 0; i <=divs; ++i)
    {
        double t = low + (i/double(divs))*(high-low);
       // double sum1 = 0.0;
        double sum2 = 0.0;

        for (int j=1; j<nn+1; ++j){
            for (double kk=0; kk<j-1; ++kk){
                auto a =
                        [&](int r, double t){
                          return r*kernel(k, t, theta, c, sig, base);
                        };
           //     sum1 += beta(n, m, a, t)* gamma_n(t-m,n);
                sum2 += beta(j, kk, a, t)*gamma_n(t-kk, j);
            }
        }
        cerr << " t= " << t << ", gamma_n(t, n)=" << gamma_n(t, n)  <<" , h_n=" << sum2 << "\n";
    }

    return 0;
}

/*
Enter kernel: PL
 Enter the value for t: 1.2
 Enter the value for n (n>=0): 6
 Enter the value for N (N>=1): 6
 Enter the value for r: 5
 Enter the value for k: 0.8
 Enter the value for theta: 0.9
 Enter the value for c: 1.1
 Enter the value for delta: 0.01
           N          n         r              rect            gamma_n              h          h_n           error

           6           2    0.000000     1.000000          0.800000          0.128942          0.000000          0.128942
           6           3    0.000000     1.000000          0.660000          0.235318          0.000000          0.235318
           6           3    1.000000     1.000000          0.660000          0.238542          0.000000          0.238542
           6           4    0.000000     1.000000          0.282667          0.284101          0.000000          0.284101
           6           4    1.000000     1.000000          0.282667          0.284316          0.000000          0.284316
           6           4    2.000000     1.000000          0.282667          0.284316          0.000000          0.284316
           6           5    0.000000     1.000000          0.086067          0.298188         -0.000000          0.298188
           6           5    1.000000     1.000000          0.086067          0.298199         -0.000000          0.298199
           6           5    2.000000     1.000000          0.086067          0.298199         -0.000000          0.298199
           6           5    3.000000     1.000000          0.086067          0.298199         -0.000000          0.298199
           6           6    0.000000     1.000000          0.020720          0.301538          0.000000          0.301538
           6           6    1.000000     1.000000          0.020720          0.301539          0.000000          0.301539
           6           6    2.000000     1.000000          0.020720          0.301539          0.000000          0.301539
           6           6    3.000000     1.000000          0.020720          0.301539          0.000000          0.301539
           6           6    4.000000     1.000000          0.020720          0.301539          0.000000          0.301539
 */