#include <Matrix.h>
#include <MatrixOperations.h>

#include <iostream>
using namespace std;



int main()
{
    int M = 7;
    int N = 4;

    Matrix2d A(M, N);
    Matrix2d B(M, N);

    for(int i = 0; i < M * N; i++)
    {
        A.Data()[i] = (double)i;
    }
    A.print();
    (+A).print();


    B = A;
    B = (+A);
    B.print();
    return 0;
}