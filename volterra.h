//
// Created by leanne on 10/25/20.
//

#ifndef EXPLICIT_VOLTERRA_VOLTERRA_H
#define EXPLICIT_VOLTERRA_VOLTERRA_H

#include <string>
#include <functional>

namespace volterra {
    double rect(double t);

    double findError(double h1,double h2);
    double gamma_n(double t, int n);
    //double beta(int n, int r, std::function<double(int r)> a);
    //void norms(double* x, int n, double& norm1, double& norminf);
    template <typename A>
    static double ar(int r, double /*t*/, const A &a)
    {
        return a(r);
    }
    template<class Integer, class Function, class Double>
    double beta(Integer n, int r, Function &&a, Double t)
    {
        if (n==0 || r==0){
            return 1;
        }
        else if (n==1 || r>=0){
            return ar(r, t, a);
        }
        else {
            double value = 0;
            for (auto l = 0; l <= n; l++)
            {
                value +=
                        beta(n, r, a, t) +
                        ar(l, t, a)*beta(n-l, r, a, t);
            }
            return value;
        }
    }
}

#endif //EXPLICIT_VOLTERRA_VOLTERRA_H
