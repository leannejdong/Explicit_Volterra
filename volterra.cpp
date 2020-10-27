//
// Created by leanne on 10/25/20.
//
#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cassert>
#include <ranges>
#include <numeric>
using namespace std;
using namespace std::ranges;
using std::cerr;

// Returns factorial of n
static int fact(int n)
{
    std::ranges::iota_view v{1, n+1};
    return std::accumulate(v.begin(), v.end(), 1, std::multiplies<>());
}


// Function definition
static int nCr(int n, int r)
{
    return fact(n) / (fact(r) * fact(n - r));
}

//static double my_fct(double t, double a, double k){
//    return (t>= a ? pow(t-a, k) : 0);
//}

static int step(double x)
{
    if (x<0) return 0;
    return 1;
}

template <typename A>
static double ar(int r, double t, const A &a)
{
    return a(r, t);
}

double
volterra::beta(
        int n,
        int r,
        std::function<double(int r, double t)> a,
        double t
)
{
    double value = 0;
    int l;

    if (n=0, r=0){
        return 1;
    }
    else if (n=0, r>=1){
        return 0;
    }
    else if (n=1, r>=0){
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

double volterra::rect(double t, double low, double up){
    return ((t>=low)  && (t< up)? 1 : 0);
}

void volterra::norms(double* x, int n, double& norm1, double& norminf)
{
    std::transform(x, x+n, x, [](double x){ return fabs(x); });
    norm1 = std::accumulate(x, x+n, 0.0);
    norminf = *std::max_element(x, x+n);
}

double volterra::findError(double h1,double h2)
{
    return fabs(h1-h2);
}

double volterra::gamma_n(double t, int n)
{
    using std::pow;
    using std::max;
    double sum = 0;

    for (int r = 0; r <= n; ++r)
    {
        double term = pow(-1, r)*nCr(n, r)*(1.0/fact(n-1))*pow(max(t-r, 0.0), n-1);
        sum += term;
    }
    return sum;
}