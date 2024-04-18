#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
using namespace std;


int main()
{

    double x;
    x = MathematicalOperations::getInstance()->log10(156);

    printf("x = %5.14lf\n", x);


    return 0;
}