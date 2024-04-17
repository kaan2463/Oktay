#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
using namespace std;


int main()
{

    double x;
    x = MathematicalOperations::getInstance()->arctan(0.707106);

    printf("x = %5.14lf\n", x);


    return 0;
}