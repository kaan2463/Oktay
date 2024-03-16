#include <MatrixOperations.h>

#include <iostream>
using namespace std;



int main()
{
    int M = 3;
    int N = 6;

    double* A = new double[M * N];
    double* B = new double[M * N];

    for(int i = 0; i < M * N; i++)
    {
        A[i] = (double)i;
    }

    MatrixOperations::getInstance()->print2d(A, M, N);

    MatrixOperations::getInstance()->transpose2d(A, M, M);

    MatrixOperations::getInstance()->print2d(A, N, M);

    return 0;
}