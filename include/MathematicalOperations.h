#pragma once
#ifndef OKTAY_MATHEMATICAL_OPERATIONS
#define OKTAY_MATHEMATICAL_OPERATIONS

#define DEFAULT_TYPE STD

typedef double (*DFUNC1D)(double);
typedef double (*DFUNC2D)(double, double);

enum MathImplementaionType
{
    BASE, STD
};


class MathematicalOperations
{
private:
    static MathematicalOperations* INSTANCE;
protected:
    MathematicalOperations() {}
public:
    static MathematicalOperations* getInstance(MathImplementaionType type = DEFAULT_TYPE);
    static MathematicalOperations* getInstanceBase();
    ~MathematicalOperations() {}

    /*
    * Absolute value
    * return |x|
    */
    virtual double abs(double x);

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
    virtual double sqrt(double x);

    /*
    * Numerical Derivative of y = f(x)
    */
    double derivative(DFUNC1D f, double x);

    /*
    * Power 10
    * return x^10
    */
    virtual double pow10(double x);

    /*
    * Power of inverse 10
    * by Newton Raphson Algorithm
    * return x^(1/10)
    */
    virtual double powInv10(double x);


    /*
    * Power
    * n : double
    * return x^n
    */
    virtual double pow(double x, double n);

    /*
    * Exponential function
    * n : double
    * return e^x
    */
    virtual double exp(double x);

    /*
    * Sinus Function in radyan
    * by using Taylor (or McLauren) Series
    * sin(x)=sum _{n=0}^{infty}{frac{(-1)^{n}}{(2n+1)!}}x^{2n+1}
    * All x
    */
    virtual double sin(double x);

    /*
    * Cosinus Function in radyan
    * by using Taylor (or McLauren) Series
    * cos(x)=sum _{n=0}^{infty }{frac {(-1)^{n}}{(2n)!}}x^{2n}
    * All x
    */
    virtual double cos(double x);

    /*
    * Tangent Function in radyan
    * by using Taylor (or McLauren) Series
    * tan(x)=sum _{n=1}^{infty }{frac {B_{2n}(-4)^{n}(1-4^{n})}{(2n)!}}x^{2n-1}
    * abs(x) <= PI / 2
    */
    virtual double tan(double x);

    /*
    * Arcsin Function in radyan
    * by using Taylor (or McLauren) Series
    * arcsin(x)=sum _{n=0}^{infty }{frac {(2n)!}{4^{n}(n!)^{2}(2n+1)}}x^{2n+1}
    * abs(x) <= 1
    */
    virtual double arcsin(double x);

    /*
    * (PI / 2) - arcsin(x)
    */
    virtual double arccos(double x);

    /*
    * arcangent Function in radyan
    * by using Taylor (or McLauren) Series
    * arctan(x)=sum _{n=0}^{infty }{frac {(-1)^{n}}{2n+1}}x^{2n+1}
    * abs(x) <= 1
    */
    virtual double arctan(double x);

    /*
    * Natural Algorithm byusing Taylor Series
    * ln(x) = sum _{n=1}^{infty }{frac {(-1)^{n+1}{z-1}^n}{n}}
    * 0 < x <= 2   >>> x = (0,inf]
    */
    virtual double ln(double x);

    /*
    * Logarithm
    * c = log_{base}(x)
    * log_{base}(x) = ln(x) / ln(base)
    */
    virtual double log(double base, double x);

    /*
    * Logarithm on base2
    * log2(x) = ln(x) / ln(2)
    */
    virtual double log2(double x);

    /*
    * Logarithm on base10
    * log10(x) = ln(x) / ln(10)
    */
    virtual double log10(double x);

};

#endif