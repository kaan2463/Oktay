#include <Matrix.h>
#include <MatrixOperations.h>

#include <iostream>
using namespace std;



int main()
{
    int M = 5;

    Matrix2d A(M, M);
    Matrix2d B(M, M);

    for(size_t i = 0; i < M * M; i++)
    {
        A.Data()[i] = (double)i + 1;

        if(i / M == i % M)
        {
            B.Data()[i] = 1;
        }
        else
        {
            B.Data()[i] = 0;
        }
    }
    A.print();
    B.print();

    MatrixOperations::getInstance()->inverse(A.Data(), B.Data(), M);

    A.print();
    B.print();

    Matrix2d C(M, M);

    C = A * B;
    C.print();

    return 0;
}