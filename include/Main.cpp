#include <Matrix.h>
#include <MatrixOperations.h>

#include <iostream>
using namespace std;



int main()
{
    int M = 2;

    // double MAT[] = {
    //     1,    2,   45,    2,
    //     1,    3,    5,    1,
    //     5,    2,    6,    2,
    //     4,    7,    3,    2 };

    double MAT[] = {
        0, 2,
        0, 4
    };

    Matrix2d A(M, M);
    Matrix2d B(M, M);

    for(size_t i = 0; i < M * M; i++)
    {
        //A.Data()[i] = (double)i * i + 1;

        A.Data()[i] = MAT[i];

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