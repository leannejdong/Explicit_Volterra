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
                cerr << ", beta_n(n, r, a, t)" << beta(j, r, a, t) << "\n";
            }
        }
        cerr << " t= " << t  << ", gamma_n(t, n)=" << gamma_n(t, n)  <<" , h_n=" << sum2 << "\n";
    }

    return 0;
}
/*
 * gamma_n(1200)=0
 t= 0, gamma_n(t, n)=0 , h_n=0
 t= 0.8, gamma_n(t, n)=2.46872e-125 , h_n=1.2224
 t= 1.6, gamma_n(t, n)=1.49225e-101 , h_n=1.92473
 t= 2.4, gamma_n(t, n)=1.21633e-87 , h_n=1.3114
 t= 3.2, gamma_n(t, n)=9.02008e-78 , h_n=0.4374
 t= 4, gamma_n(t, n)=4.08415e-70 , h_n=0.0416667
 t= 4.8, gamma_n(t, n)=7.35224e-64 , h_n=6.66667e-05
 t= 5.6, gamma_n(t, n)=1.42959e-58 , h_n=0
 t= 6.4, gamma_n(t, n)=5.45166e-54 , h_n=0
 t= 7.2, gamma_n(t, n)=5.98924e-50 , h_n=0
 t= 8, gamma_n(t, n)=2.46354e-46 , h_n=0
 t= 8.8, gamma_n(t, n)=4.57049e-43 , h_n=0
 t= 9.6, gamma_n(t, n)=4.38448e-40 , h_n=0
 t= 10.4, gamma_n(t, n)=2.41038e-37 , h_n=0
 t= 11.2, gamma_n(t, n)=8.21886e-35 , h_n=0
 t= 12, gamma_n(t, n)=1.84923e-32 , h_n=0
 t= 12.8, gamma_n(t, n)=2.88426e-30 , h_n=0
 t= 13.6, gamma_n(t, n)=3.24512e-28 , h_n=0
 t= 14.4, gamma_n(t, n)=2.72098e-26 , h_n=0
 t= 15.2, gamma_n(t, n)=1.74675e-24 , h_n=0
 t= 16, gamma_n(t, n)=8.78097e-23 , h_n=0
 t= 16.8, gamma_n(t, n)=3.52323e-21 , h_n=0
 t= 17.6, gamma_n(t, n)=1.14678e-19 , h_n=0
 t= 18.4, gamma_n(t, n)=3.07053e-18 , h_n=0
 t= 19.2, gamma_n(t, n)=6.84516e-17 , h_n=0
 t= 20, gamma_n(t, n)=1.28395e-15 , h_n=0
 t= 20.8, gamma_n(t, n)=2.045e-14 , h_n=0
 t= 21.6, gamma_n(t, n)=2.78826e-13 , h_n=0
 t= 22.4, gamma_n(t, n)=3.27765e-12 , h_n=0
 t= 23.2, gamma_n(t, n)=3.34297e-11 , h_n=0
 t= 24, gamma_n(t, n)=2.97499e-10 , h_n=0
 t= 24.8, gamma_n(t, n)=2.32168e-09 , h_n=0
 t= 25.6, gamma_n(t, n)=1.59601e-08 , h_n=0
 t= 26.4, gamma_n(t, n)=9.70349e-08 , h_n=0
 t= 27.2, gamma_n(t, n)=5.23659e-07 , h_n=0
 t= 28, gamma_n(t, n)=2.51653e-06 , h_n=0
 t= 28.8, gamma_n(t, n)=1.08007e-05 , h_n=0
 t= 29.6, gamma_n(t, n)=4.15072e-05 , h_n=0
 t= 30.4, gamma_n(t, n)=0.000143164 , h_n=0
 t= 31.2, gamma_n(t, n)=0.00044411 , h_n=0
 t= 32, gamma_n(t, n)=0.00124136 , h_n=0
 t= 32.8, gamma_n(t, n)=0.0031316 , h_n=0
 t= 33.6, gamma_n(t, n)=0.00714057 , h_n=0
 t= 34.4, gamma_n(t, n)=0.0147348 , h_n=0
 t= 35.2, gamma_n(t, n)=0.0275473 , h_n=0
 t= 36, gamma_n(t, n)=0.0467025 , h_n=0
 t= 36.8, gamma_n(t, n)=0.0718566 , h_n=0
 t= 37.6, gamma_n(t, n)=0.100397 , h_n=0
 t= 38.4, gamma_n(t, n)=0.127438 , h_n=0
 t= 39.2, gamma_n(t, n)=0.146967 , h_n=0
 t= 40, gamma_n(t, n)=0.154001 , h_n=0
 t= 40.8, gamma_n(t, n)=0.146967 , h_n=0
 t= 41.6, gamma_n(t, n)=0.127438 , h_n=0
 t= 42.4, gamma_n(t, n)=0.100397 , h_n=0
 t= 43.2, gamma_n(t, n)=0.0718566 , h_n=0
 t= 44, gamma_n(t, n)=0.0467025 , h_n=0
 t= 44.8, gamma_n(t, n)=0.0275473 , h_n=0
 t= 45.6, gamma_n(t, n)=0.0147348 , h_n=0
 t= 46.4, gamma_n(t, n)=0.00714057 , h_n=0
 t= 47.2, gamma_n(t, n)=0.0031316 , h_n=0
 t= 48, gamma_n(t, n)=0.00124136 , h_n=0
 t= 48.8, gamma_n(t, n)=0.00044411 , h_n=0
 t= 49.6, gamma_n(t, n)=0.000143164 , h_n=0
 t= 50.4, gamma_n(t, n)=4.15072e-05 , h_n=0
 t= 51.2, gamma_n(t, n)=1.08007e-05 , h_n=0
 t= 52, gamma_n(t, n)=2.51653e-06 , h_n=0
 t= 52.8, gamma_n(t, n)=5.23659e-07 , h_n=0
 t= 53.6, gamma_n(t, n)=9.70349e-08 , h_n=0
 t= 54.4, gamma_n(t, n)=1.59601e-08 , h_n=0
 t= 55.2, gamma_n(t, n)=2.32168e-09 , h_n=0
 t= 56, gamma_n(t, n)=2.97499e-10 , h_n=0
 t= 56.8, gamma_n(t, n)=3.34297e-11 , h_n=0
 t= 57.6, gamma_n(t, n)=3.27765e-12 , h_n=0
 t= 58.4, gamma_n(t, n)=2.78826e-13 , h_n=0
 t= 59.2, gamma_n(t, n)=2.045e-14 , h_n=0
 t= 60, gamma_n(t, n)=1.28395e-15 , h_n=0
 t= 60.8, gamma_n(t, n)=6.84516e-17 , h_n=0
 t= 61.6, gamma_n(t, n)=3.07053e-18 , h_n=0
 t= 62.4, gamma_n(t, n)=1.14678e-19 , h_n=0
 t= 63.2, gamma_n(t, n)=3.52323e-21 , h_n=0
 t= 64, gamma_n(t, n)=8.78097e-23 , h_n=0
 t= 64.8, gamma_n(t, n)=1.74675e-24 , h_n=0
 t= 65.6, gamma_n(t, n)=2.72098e-26 , h_n=0
 t= 66.4, gamma_n(t, n)=3.24512e-28 , h_n=0
 t= 67.2, gamma_n(t, n)=2.88426e-30 , h_n=0
 t= 68, gamma_n(t, n)=1.84923e-32 , h_n=0
 t= 68.8, gamma_n(t, n)=8.21886e-35 , h_n=0
 t= 69.6, gamma_n(t, n)=2.41038e-37 , h_n=0
 t= 70.4, gamma_n(t, n)=4.38448e-40 , h_n=0
 t= 71.2, gamma_n(t, n)=4.57049e-43 , h_n=0
 t= 72, gamma_n(t, n)=2.46354e-46 , h_n=0
 t= 72.8, gamma_n(t, n)=5.98924e-50 , h_n=0
 t= 73.6, gamma_n(t, n)=5.45166e-54 , h_n=0
 t= 74.4, gamma_n(t, n)=1.42959e-58 , h_n=0
 t= 75.2, gamma_n(t, n)=7.35224e-64 , h_n=0
 t= 76, gamma_n(t, n)=4.08415e-70 , h_n=0
 t= 76.8, gamma_n(t, n)=9.02008e-78 , h_n=0
 t= 77.6, gamma_n(t, n)=1.21633e-87 , h_n=0
 t= 78.4, gamma_n(t, n)=1.49225e-101 , h_n=0
 t= 79.2, gamma_n(t, n)=2.46872e-125 , h_n=0
 t= 80, gamma_n(t, n)=0 , h_n=0
 */
