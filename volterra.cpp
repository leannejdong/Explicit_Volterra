//
// Created by leanne on 10/25/20.
//

#include "volterra.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cassert>
using namespace std;


// Returns factorial of n
static int fact(int n)
{
    int res = 1;
    for (int i = 2; i <= n; i++)
        res = res * i;
    return res;
}


// Function definition
static int nCr(int n, int r)
{
    return fact(n) / (fact(r) * fact(n - r));
}

static double
rect(double t, double low, double up){
    return ((t>=low)  && (t< up)? 1 : 0);
}


static double impul(double t){
    return (t==0 ? 1 : 0);
}


/*double taPow(double a, double r){
    return pow(t-a, r)
}*/


static double my_fct(double t, double a, double k){
    return (t>= a ? pow(t-a, k) : 0);
}


static double
kernel(
        int r,
        double k,
        double t,
        double theta,
        double c,
        double sig,
        double del,
        std:: string& base
)
{
    if(base== "PL"){
        return k*theta*pow(c,theta)/(pow((c+(r+1)*del),(1+theta)));
    }
    else if (base=="Exp"){
        return k*theta*exp(-theta*t);
    }
    else if (base=="RL"){
        return (k*t/(pow(sig,2)))*exp(-t*t/(2*pow(sig,2)));
    }

    assert(false); // shouldn't happen
}


static double
ar(
        int r,
        double k,
        double t,
        double theta,
        double c,
        double sig,
        double del,
        std:: string& base
)
{
    return del*kernel(r,  k,  t,  theta,  c,  sig,  del, base);
}


double
volterra::beta(
        int n,
        int r,
        double k,
        double t,
        double theta,
        double c,
        double sig,
        double del,
        std:: string& base
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
        return ar(r,k,t,theta,c,sig,del,base);
    }
    else {
        value =
                beta(n,r,k,t,theta,c,sig,del,base) +
                ar(1,k,t,theta,c,sig,del,base)*beta(n-1,r,k,t,theta,c,sig,del,base);

        // lacking an argument, also need to incorporate kernel
        return value;
    }
}


void volterra::norms(double* x, int n, double& norm1, double& norminf)
{
    norm1 = fabs(x[0]);
    norminf = norm1;
    for(int i=1; i<n; i++)
    {
        x[i] = fabs(x[i]);
        norm1 += x[i];
        if(x[i]> norminf)
            norminf = x[i];
    }
}


double volterra::findError(double h1,double h2)
{
    return fabs(h1-h2);
}


double volterra::series(int n, double t){
// calculating the value of (n-1)!

    int n1fact = fact(n-1);
    double sum = 0.0;

// loop to display the series
    for ( int i=0; i<n+1; i++) {

        int mPow = pow(-1,i);
        int comb = nCr(n, i);
        double myfct = my_fct(t,i,n-1);
        sum += mPow * comb * myfct/ n1fact;
    }
    return sum;
}