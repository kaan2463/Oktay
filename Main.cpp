#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
using namespace std;

int main()
{

    double x;
    x = MathematicalOperations::getInstance()->pow(0.49, -0.5);

    printf("x = %5.14lf\n", x);


    return 0;
}