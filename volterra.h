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
    double beta(int n, int r, std::function<double(int r)> a);
    void norms(double* x, int n, double& norm1, double& norminf);
    double beta(int n, int r, std::function<double(int)> a, double t);
}

#endif //EXPLICIT_VOLTERRA_VOLTERRA_H
