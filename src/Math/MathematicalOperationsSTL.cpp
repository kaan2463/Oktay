#include <cmath>
#include <MathematicalOperationsSTL.h>

#include <Exception.h>

typedef unsigned long long uint64_t;

MathematicalOperationsSTL* MathematicalOperationsSTL::INSTANCE = 0;

MathematicalOperations* MathematicalOperationsSTL::getInstance()
{
    if(INSTANCE == 0)
    {
        INSTANCE = new MathematicalOperationsSTL;
    }
    return INSTANCE;
}

double MathematicalOperationsSTL::abs(double x)
{
    return std::abs(x);
}

double MathematicalOperationsSTL::sqrt(double x)
{
    return std::sqrt(x);
}

double MathematicalOperationsSTL::pow10(double x)
{
    return std::pow(x, 10);
}

double MathematicalOperationsSTL::powInv10(double x)
{
    return std::pow(x, 0.1);
}

double MathematicalOperationsSTL::pow(double x, double n)
{
    return std::pow(x, n);
}

double MathematicalOperationsSTL::exp(double x)
{
    return std::exp(x);
}

double MathematicalOperationsSTL::sin(double x)
{
    return std::sin(x);
}

double MathematicalOperationsSTL::cos(double x)
{
    return std::cos(x);
}


double MathematicalOperationsSTL::tan(double x)
{
    return std::tan(x);
}

double MathematicalOperationsSTL::arcsin(double x)
{
    return std::asin(x);
}

double MathematicalOperationsSTL::arccos(double x)
{
    return std::acos(x);
}

double MathematicalOperationsSTL::arctan(double x)
{
    return std::atan(x);
}


double MathematicalOperationsSTL::ln(double x)
{
    return std::log(x);
}

double MathematicalOperationsSTL::log(double base, double x)
{
    return std::log(base) / std::log(x);
}

double MathematicalOperationsSTL::log2(double x)
{
    return std::log2(x);
}

double MathematicalOperationsSTL::log10(double x)
{
    return std::log10(x);
}