//
// Created by leanne on 10/25/20.
//
#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cassert>
//#include <ranges>
#include <numeric>
using namespace std;
//using namespace std::ranges;
using std::cerr;

// Returns factorial of n
//static double fact(int n)
//{
//    std::ranges::iota_view v{1, n+1};
//    return std::accumulate(v.begin(), v.end(), 1, std::multiplies<>());
//}

static double fact(int n)
{
    double product = 1;

    while (n > 1) {
        product *= n;
        --n;
    }

    return product;
}

template <typename A>
static double ar(int r, double /*t*/, const A &a)
{
    return a(r);
}

double
volterra::beta(int n, int r, std::function<double(int r)> a, double t)
{
    double value = 0;
    int l;

    if (n==0 || r==0){
        return 1;
    }
    else if (n==1 || r>=0){
        return ar(r, t, a);
    }
    else {
        for (l = 0; l <= n; l++)
        {
            value +=
                    beta(n, r, a, t) +
                    ar(l, t, a)*beta(n-l, r, a, t);
        }
        return value;
    }
}
#if 0
double volterra::rect(double t){
    return ((t>=1)  && (t< 0)? 1 : 0);
}
#endif

//void volterra::norms(double* x, int n, double& norm1, double& norminf)
//{
//    std::transform(x, x+n, x, [](double x){ return fabs(x); });
//    norm1 = std::accumulate(x, x+n, 0.0);
//    norminf = *std::max_element(x, x+n);
//}

double volterra::findError(double h1,double h2)
{
    return fabs(h1-h2);
}

double volterra::gamma_n(double t, int n)
{
    if (t > n/2) {
        t = n-t;
    }
    using std::pow;
    using std::max;
    double sum = 0;

    for (int r = 0; r <= n; ++r)
    {
       // double term = pow(-1, r)*nCr(n, r)*(1.0/fact(n-1))*pow(max(t-r, 0.0), n-1);
        double term = pow(-1,r)*(n/(fact(r)*fact(n-r)))*pow(max(t-r,0.0),n-1);
        sum += term;
    }
    return sum;
}