#include <MathematicalOperations.h>

typedef unsigned long long uint64_t;

MathematicalOperations* MathematicalOperations::INSTANCE = 0;

MathematicalOperations* MathematicalOperations::getInstance()
{
    if(INSTANCE == 0)
    {
        INSTANCE = new MathematicalOperations;
    }
    return INSTANCE;
}

inline double abs(double x)
{
    // breaks strict aliasing, but compiler writer knows this behavior for the platform
    uint64_t i = reinterpret_cast<const uint64_t&>(x);
    i &= 0x7FFFFFFFFFFFFFFFULL; // clear sign bit

    return reinterpret_cast<const double&>(i);
}


void MathematicalOperations::newtonRaphson(DFUNC2D f, DFUNC2D d, double initialPoint, double* result)
{
    size_t iter = 0;
    *result = initialPoint;
    double old_result;
    do
    {
        old_result = *result;
        *result = *result - (f(*result) / d(*result));
        iter++;
    } while(iter < MAX_ITER && abs(old_result - (*result)) > DE);
}

void MathematicalOperations::newtonRaphson(DFUNC2DA f, DFUNC2DA d, double initialPoint, double* result, double value)
{
    size_t iter = 0;
    *result = initialPoint;
    double old_result;
    do
    {
        old_result = *result;
        *result = *result - (f(*result, value) / d(*result, value));
        iter++;
    } while(iter < MAX_ITER && abs(old_result - (*result)) > DE);
}

inline double parabola(double x, double a)
{
    return x * x - a;
}

inline double derivParabola(double x, double a)
{
    return 2 * x;
}

double MathematicalOperations::sqrt(double x)
{
    double result;
    newtonRaphson(parabola, derivParabola, 1.0, &result, x);
    return result;
}
