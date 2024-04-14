#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
using namespace std;

double sqr(double x)
{
    return x * x - 9;
}

double dsqr(double x)
{
    return 2 * x;
}

int main()
{
    double x = MathematicalOperations::getInstance()->sqrt(225);

    printf("x = %lf\n", x);


    return 0;
}