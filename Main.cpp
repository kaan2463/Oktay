#include <MathematicalOperations.h>
#include <MatrixOperations.h>

#define MAT MathematicalOperations::getInstance(MathImplementaionType::STD)


#include <iostream>
using namespace std;

int main()
{

    printf("sin = %lf\n", MAT->log2(32));

}