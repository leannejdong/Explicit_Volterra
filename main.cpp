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
        double kappa,
        double t,
        double theta,
        double c,
        double sig,
       // double,
        std:: string& base
)
{
    if(base== "PL"){
        // return kappa*theta*pow(c,theta)/(pow((c+(r+1)*del),(1+theta)));
        return kappa*theta*pow(c,theta)/(pow((c+t),(1+theta)));
    }
    else if (base=="Exp"){
        return kappa*theta*exp(-theta*t);
    }
    else if (base=="RL"){
        return (kappa*t/(pow(sig,2)))*exp(-t*t/(2*pow(sig,2)));
    }
    else {
        assert(false); // shouldn't happen
        return 0;
    }
}


int main() {
    using volterra::rect;
    using volterra::beta;
    using volterra::gamma_n;
    using volterra::findError;
//    using volterra::norms;

    int n, nn;
    double kappa, theta, c, sig;
    std::string base;

    base = "PL"; n = 80; nn = 5; kappa = 0.90; theta = 1.33; c = 1.1;

    cout << fixed << setprecision(6);

    int divs = 100;
    double low = 0;
    double high = n;
    cerr << "gamma_n(1200)=" << gamma_n(1200,n) << "\n";
    for (int i = 0; i <=divs; ++i)
    {
        double t = low + (i/double(divs))*(high-low);
        double sum2 = 0.0;

        for (int j=1; j<nn+1; ++j){
            for (double r=0; r<j-1; ++r){
                auto a =
                        [&](int r /*t*/){
                          double t = low + (r/double(divs))*(high-low);
                          return kernel(kappa, t, theta, c, sig, base);
                        };
                sum2 += beta(j, r, a, t)*gamma_n(t-r, j);
                cerr << " t= " << t << ",  beta_n(j, r, a, t)=" << beta(j, r, a, t) << ", gamma_n(t-r, j)=" << gamma_n(t-r, j)
                     << " h_n= " << sum2 << "\n";
            }
        }
        //cerr << " t= " << t  << ", gamma_n(t, n)=" << gamma_n(t, n)  <<" , h_n=" << sum2 << "\n";
    }

    return 0;
}
