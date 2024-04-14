#pragma once
#ifndef OKTAY_MATHEMATICAL_OPERATIONS
#define OKTAY_MATHEMATICAL_OPERATIONS

#define MAX_ITER    30
#define DE          1.0e-8

typedef double (*DFUNC2D)(double);
typedef double (*DFUNC2DA)(double, double);

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
    void newtonRaphson(DFUNC2D f, DFUNC2D d, double initialPoint, double* result);

    /*
    * Newton Raphson method
    */
    void newtonRaphson(DFUNC2DA f, DFUNC2DA d, double initialPoint, double* result, double value);

    /*
    * Square Root by using Newton Raphson method
    * return sqrt(x)
    */
    double sqrt(double x);
};
#endif