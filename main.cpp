#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>

using namespace std;

enum class base { PL, Exp, RL};
auto
kernel(
        //int,
        double kappa,
        double t,
        double theta,
        double c,
        double sig,
       // double,
        base foo
)
{
    switch (foo)
    {
        case base::PL:
            return kappa*theta*pow(c,theta)/(pow((c+t), (1+theta)));
        case base::Exp:
            return kappa*theta*exp(-theta*t);
        case base::RL:
            return (kappa*t/(pow(sig,2)))*exp(-t*t/(2*pow(sig,2)));

        default:
            return 0.0;
    }
}

int main() {
    using volterra::rect;
    using volterra::beta;
    using volterra::gamma_n;
//    using volterra::findError;
//    using volterra::norms;

    int n, nn;
    double kappa, theta, c, sig;
   // std::string base;

    n = 80; nn = 5;
    //kappa = 0.80; theta = 0.8; c = 10; // PL
    kappa = 0.60; theta = 0.8; // Exp

    cout << fixed << setprecision(6);

    int divs = 100;
    double low = 0;
    double high = n;
   // cerr << "gamma_n(1200)=" << gamma_n(1200,n) << "\n";
   // std::ofstream of_beta("/home/leanne/CLionProjects/Explicit-Volterra/outputs/beta_n.txt");
    std::ofstream of_beta("outputs/beta_n_exp.txt");
   // std::ofstream of_gamma("/home/leanne/CLionProjects/Explicit-Volterra/outputs/gamma_n.txt");
    std::ofstream of_gamma("outputs/gamma_n_exp.txt");
    std::ofstream of_sum2("outputs/sum2_exp.txt");
    std::ofstream of_sum1("outputs/sum1_exp.txt");
    std::ofstream of_all("outputs/all_exp.txt");

    for (int i = 0; i <=divs; ++i)
    {
        double t = low + (i/double(divs))*(high-low);
        double sum1 = 0.0;
        double sum2 = 0.0;
        sum1 += kernel(kappa, t, theta, c, sig, base::Exp);
        of_sum1 << setw(10) << t << setw(20) << sum1 << "\n";
        for (int j=1; j<nn+1; ++j){
            for (double r=0; r<j-1; ++r) {
                auto a =
                        [&](int r /*t*/) {
                            double t = low + (r / double(divs)) * (high - low);
                            return kernel(kappa, t, theta, c, sig, base::Exp);
                        };
                sum2 += beta(j, r, a, t) * gamma_n(t - r, j);
                //                cerr << " t= " << t << ",  beta_n(j,r,a,t)=" << beta(j, r, a, t) << ", gamma_n(t-r,j)=" << gamma_n(t-r, j)
                //                     << " h_n= " << sum2 << "\n";
                of_gamma << gamma_n(t-r, j) << "\n";
                of_beta << beta(j, r, a, t) << "\n";
                of_sum2 << setw(10) << t << setw(20) << sum2 << "\n";
                of_all << setw(10) << t << setw(20)<< gamma_n(t-r, j) << "," << setw(20) << beta(j, r, a, t) << setw(20) << sum2 << "\n";
            }
        }
        //cerr << " t= " << t  << ", gamma_n(t, n)=" << gamma_n(t, n)  <<" , h_n=" << sum2 << "\n";
    }

    return 0;
}
