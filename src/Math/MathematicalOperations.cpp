#include <Exception.h>
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
    uint64_t i = reinterpret_cast<const uint64_t&>(x);
    i &= 0x7FFFFFFFFFFFFFFFULL; // clear sign bit

    return reinterpret_cast<const double&>(i);
}

double MathematicalOperations::abs(double x)
{
    return ::abs(x);
}


void MathematicalOperations::newtonRaphson(DFUNC1D f, DFUNC1D d, double initialPoint, double* result)
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

void MathematicalOperations::newtonRaphson(DFUNC2D f, DFUNC2D d, double initialPoint, double* result, double value)
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

double MathematicalOperations::derivative(DFUNC1D f, double x)
{
    double h = 1.0e-5;
    h = abs(x) * h + h;
    return (f(x + h) - f(x - h)) / (2 * h);
}

double MathematicalOperations::pow10(double x)
{
    double result = x * x;
    result = result * result * x;
    return result * result;
}

inline double pow10H(double x)
{
    double result = x * x;
    result = result * result * x;
    return result * result;
}

double MathematicalOperations::powInv10(double x)
{
    size_t iter = 0;
    double result = 1.0;
    double old_result;
    do
    {
        old_result = result;
        result = result - ((pow10(result) - x) / derivative(pow10H, result));
        iter++;
    } while(iter < MAX_ITER && abs(old_result - (result)) > DE);
    return result;
}

double powH(double x, long n)
{
    if(n == 0)
    {
        return 1.0;
    }

    if(x < 0)
    {
        if(n % 2 == 1)
        {
            return -1.0 * powH(abs(x), n);
        }
        else
        {
            throw RuntimeException("NOT REAL!!!"); //NAN (not reel)
        }
    }

    if(n < 0)
    {
        return powH(1.0 / x, abs((double)n));
    }

    if(n % 2 == 1)
    {
        return x * powH(x, n - 1);
    }

    double result = powH(x, n / 2);

    return result * result;

}

double MathematicalOperations::pow(double x, double n)
{

    if(n < 0)
    {
        return pow(1.0 / x, abs(n));
    }

    double floatingPow = n - (long long)n;

    if(floatingPow == 0.0 || floatingPow < DE)
    {
        return powH(x, ((long)n));
    }
    else
    {
        if(x < 0)
        {
            throw RuntimeException("NOT REAL!!!"); //NAN (not reel)
        }
    }

    double base = x;
    double result = 1.0;
    long fPart;
    for(size_t i = 0; i < POW_PRECISION; i++)
    {
        base = powInv10(base);
        fPart = (long)(floatingPow * 10.0);
        floatingPow = floatingPow * 10.0 - (double)fPart;
        result = result * powH(base, fPart);
    }
    return result * powH(x, ((long)n));
}

double MathematicalOperations::exp(double x)
{
    return pow(OKTAY_E, x);
}

double MathematicalOperations::sin(double x)
{
    double powX = 1.0;
    double fac = 1.0;
    double result = 0.0;
    double sign;
    for(size_t i = 0; i < TRIGONOMETRIC_DEPTH - 1; i++)
    {
        if(i % 2 == 1)
        {
            sign = ((i - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
            result += sign * powX / fac;
        }
        fac *= (double)(i + 1);
        powX = powX * x;
    }
    return result;
}

double MathematicalOperations::cos(double x)
{
    double powX = 1.0;
    double fac = 1.0;
    double result = 0.0;
    double sign;
    for(size_t i = 0; i < TRIGONOMETRIC_DEPTH; i++)
    {
        if(i % 2 == 0)
        {
            sign = (i / 2) % 2 == 0 ? 1.0 : -1.0;
            result += sign * powX / fac;
        }
        fac *= (double)(i + 1);
        powX = powX * x;
    }
    return result;
}


double MathematicalOperations::tan(double x)
{
    double powX = 1.0;
    double fac = 1.0;
    double resultSin = 0.0;
    double resultCos = 0.0;
    double signSin;
    double signCos;
    for(size_t i = 0; i < TRIGONOMETRIC_DEPTH; i++)
    {
        if(i % 2 == 1)
        {
            signSin = ((i - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
            resultSin += signSin * powX / fac;
        }

        if(i % 2 == 0)
        {
            signCos = (i / 2) % 2 == 0 ? 1.0 : -1.0;
            resultCos += signCos * powX / fac;
        }
        fac *= (double)(i + 1);
        powX = powX * x;
    }
    return resultSin / resultCos;
}

double MathematicalOperations::arcsin(double x)
{
    double powX = 1.0;
    double pow4 = 1.0;
    double fac = 1.0;
    double fac2n = 1.0;
    double result = 0.0;
    for(size_t i = 0; i < TRIGONOMETRIC_DEPTH + 1; i++)
    {
        if(i > 1)
        {
            fac2n *= (double)(i - 1);
        }
        if(i % 2 == 1)
        {
            if(i > 1)
            {
                fac *= ((double)(i - 1)) / 2.0;
            }
            fac = fac == 0 ? 1.0 : fac;
            result += (fac2n * powX) / (pow4 * (fac * fac) * ((double)i));
            pow4 *= 4;
        }
        powX = powX * x;
    }
    return result;
}

double MathematicalOperations::arccos(double x)
{
    return (OKTAY_PI / 2.0) - arcsin(x);
}

double MathematicalOperations::arctan(double x)
{
    double powX = 1.0;
    double result = 0.0;
    double sign;
    for(size_t i = 0; i < TRIGONOMETRIC_DEPTH - 1; i++)
    {
        if(i % 2 == 1)
        {
            sign = ((i - 1) / 2) % 2 == 0 ? 1.0 : -1.0;

            result += sign * powX / ((double)i);
        }
        powX = powX * x;
    }
    return result;
}


double MathematicalOperations::ln(double x)
{
    double result = 0.0;
    double sign;
    double z = x > 2 ? (1 / x) - 1 : x - 1;
    double powX = z;
    for(size_t i = 1; i < LOGARITMIC_DEPTH; i++)
    {
        sign = i % 2 == 1 ? 1.0 : -1.0;
        result += (sign * powX) / ((double)i);
        powX *= z;
    }
    return x > 2 ? -result : result;
}

double MathematicalOperations::log(double base, double x)
{
    double resultBase = 0.0;
    double resultX = 0.0;
    double sign;
    double zX = x > 2 ? (1 / x) - 1 : x - 1;
    double zBase = base > 2 ? (1 / base) - 1 : base - 1;
    double powBase = zBase;
    double powX = zX;
    for(size_t i = 1; i < LOGARITMIC_DEPTH; i++)
    {
        sign = i % 2 == 1 ? 1.0 : -1.0;
        resultX += (sign * powX) / ((double)i);
        resultBase += (sign * powBase) / ((double)i);
        powX *= zX;
        powBase *= zBase;
    }

    resultBase = base > 2 ? -resultBase : resultBase;
    resultX = x > 2 ? -resultX : resultX;

    return (resultX / resultBase);
}

double MathematicalOperations::log2(double x)
{
    return log(2.0, x);
}

double MathematicalOperations::log10(double x)
{
    return log(10.0, x);
}