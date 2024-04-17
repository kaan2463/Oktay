#pragma once
#ifndef OKTAY_MATHEMATICAL_OPERATIONS
#define OKTAY_MATHEMATICAL_OPERATIONS

#define MAX_ITER    150
#define DE          1.0e-8

#define POW_PRECISION 8

#define TRIGONOMETRIC_DEPTH 19

#define OKTAY_PI 3.14159265358979323846

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

    /*
    * Sinus Function in radyan
    * by using Taylor (or McLauren) Series
    * sin(x)=sum _{n=0}^{infty}{frac{(-1)^{n}}{(2n+1)!}}x^{2n+1}
    * All x
    */
    double sin(double x);

    /*
    * Cosinus Function in radyan
    * by using Taylor (or McLauren) Series
    * cos(x)=sum _{n=0}^{infty }{frac {(-1)^{n}}{(2n)!}}x^{2n}
    * All x
    */
    double cos(double x);

    /*
    * Tangent Function in radyan
    * by using Taylor (or McLauren) Series
    * tan(x)=sum _{n=1}^{infty }{frac {B_{2n}(-4)^{n}(1-4^{n})}{(2n)!}}x^{2n-1}
    * abs(x) <= PI / 2
    */
    double tan(double x);

    /*
    * Arcsin Function in radyan
    * by using Taylor (or McLauren) Series
    * arcsin(x)=sum _{n=0}^{infty }{frac {(2n)!}{4^{n}(n!)^{2}(2n+1)}}x^{2n+1}
    * abs(x) <= 1
    */
    double arcsin(double x);

    /*
    * (PI / 2) - arcsin(x)
    */
    double arccos(double x);

    /*
    * arcangent Function in radyan
    * by using Taylor (or McLauren) Series
    * arctan(x)=sum _{n=0}^{infty }{frac {(-1)^{n}}{2n+1}}x^{2n+1}
    * abs(x) <= 1
    */
    double arctan(double x);


};
#endif