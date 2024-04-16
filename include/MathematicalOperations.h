#pragma once
#ifndef OKTAY_MATHEMATICAL_OPERATIONS
#define OKTAY_MATHEMATICAL_OPERATIONS

#define MAX_ITER    150
#define DE          1.0e-8

#define POW_PRECISION 8

typedef double (*DFUNC1D)(double);
typedef double (*DFUNC2D)(double, double);

class MathematicalOperations
{
private:
    static MathematicalOperations* INSTANCE;
    MathematicalOperations() {}
public:
    static MathematicalOperations* getInstance();
    ~MathematicalOperations() {}

    /*
    * Newton Raphson method
    * x0 <- initialPoint
    * x = x0 - f/f'
    * f(result) = 0
    */
    void newtonRaphson(DFUNC1D f, DFUNC1D d, double initialPoint, double* result);

    /*
    * Newton Raphson method
    */
    void newtonRaphson(DFUNC2D f, DFUNC2D d, double initialPoint, double* result, double value);

    /*
    * Square Root by using Newton Raphson method
    * return sqrt(x)
    */
    double sqrt(double x);

    /*
    * Numerical Derivative of y = f(x)
    */
    double derivative(DFUNC1D f, double x);

    /*
    * Power 10
    * return x^10
    */
    double pow10(double x);

    /*
    * Power of inverse 10
    * by Newton Raphson Algorithm
    * return x^(1/10)
    */
    double powInv10(double x);


    /*
    * Power
    * n : double
    * return x^n
    */
    double pow(double x, double n);
};
#endif