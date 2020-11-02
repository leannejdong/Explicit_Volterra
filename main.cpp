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

double h(double t, const std::vector<double> &ts, double kappa, double theta)
{
    using std::exp;
    double sum1 = 0;
    std::ofstream of_sum1("outputs/sum1_exp.txt");
    for (double ti : ts) {
        sum1 += kappa*theta*exp(-theta*(t-ti));
        of_sum1 << setw(10) << ti << setw(20) << sum1 << "\n";
    }
    return sum1;
}

void test(double kappa, double theta)
{
    cout << fixed << setprecision(6);
    double low = 0;
    double high = 10;
    int divs = 100;
    std::vector<double> ts;

    for (int i=0; i<=divs; ++i) {
        double t = low + (i/double(divs))*(high-low);
        ts.push_back(t);
    }

    for (double t: ts) {
        h(t, ts, kappa, theta);
    }
}

int main() {
    using volterra::rect;
    using volterra::beta;
    using volterra::gamma_n;
//    using volterra::findError;
//    using volterra::norms;

    int /*n,*/ nn;
    double kappa, theta, c, sig;
   // std::string base;

    /*n = 10;*/ nn = 5;
    //kappa = 0.80; theta = 0.8; c = 10; // PL
    kappa = 0.60; theta = 0.8; // Exp

    cout << fixed << setprecision(6);
    double t = 1.0;
    int divs = 100;
    double low = 0;
    double high = 1.0;
    std::vector<double> ts;
    for (int i=0; i<=divs; ++i) {
        double ti = low + (i/double(divs))*(high-low);
        ts.push_back(ti);
    }
    h(t, ts, kappa, theta);


    std::ofstream of_beta("outputs/beta_n_exp.txt");
    std::ofstream of_gamma("outputs/gamma_n_exp.txt");
    std::ofstream of_sum2("outputs/sum2_exp.txt");
    std::ofstream of_all("outputs/all_exp.txt");

    for (int i = 0; i <=divs; ++i)
    {
        double ti = low + (i/double(divs))*(high-low);
        double sum2 = 0.0;
        for (int j=1; j<nn+1; ++j){
            for (double r=0; r<j-1; ++r) {
                auto a =
                        [&](int r /*t*/) {
                            double ti = low + (r / double(divs)) * (high - low);
                            return kernel(kappa, ti, theta, c, sig, base::Exp);
                        };
                sum2 += beta(j, r, a, ti) * gamma_n(ti - r, j);
                //                cerr << " t= " << t << ",  beta_n(j,r,a,t)=" << beta(j, r, a, t) << ", gamma_n(t-r,j)=" << gamma_n(t-r, j)
                //                     << " h_n= " << sum2 << "\n";
                of_gamma << gamma_n(ti-r, j) << "\n";
                of_beta << beta(j, r, a, ti) << "\n";
                of_sum2 << setw(10) << ti << setw(20) << sum2 << "\n";
                of_all << setw(10) << ti << setw(20)<< gamma_n(ti-r, j) << "," << setw(20) << beta(j, r, a, ti) << setw(20) << sum2 << "\n";
            }
        }
        //cerr << " t= " << t  << ", gamma_n(t, n)=" << gamma_n(t, n)  <<" , h_n=" << sum2 << "\n";
    }


    return 0;
}
