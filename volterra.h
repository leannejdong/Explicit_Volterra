//
// Created by leanne on 10/25/20.
//

#ifndef EXPLICIT_VOLTERRA_VOLTERRA_H
#define EXPLICIT_VOLTERRA_VOLTERRA_H

#include <string>

namespace volterra {
    double series(int n, double t);

    double findError(double h1,double h2);

    double
    beta(
            int n,
            int r,
            double k,
            double t,
            double theta,
            double c,
            double sig,
            double del,
            std:: string& base
    );

    void norms(double* x, int n, double& norm1, double& norminf);
};

#endif //EXPLICIT_VOLTERRA_VOLTERRA_H
