#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
#include <MatrixOperations.h>
using namespace std;


int main()
{

    size_t M = 5;
    size_t N = 4;

    //double* A = new double[M * N];

    double A[] = { 1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 5,6,7,8 };
    double B[] = {
        3,4,5,
        4,5,6,
        5,6,7,
        6,7,8 };

    double* Q = new double[M * M];
    double* R = new double[M * N];

    double* X = new double[M * N];

    //for(size_t i = 0; i < M; i++)
    //{
    //    for(size_t j = 0; j < N; j++)
    //    {
    //        A[i * N + j] = (i + 0.02 * i + 1.0);
    //    }
    //}

    MatrixOperations::getInstance()->qr(A, Q, R, M, N);

    MatrixOperations::getInstance()->print2d(A, M, N);
    MatrixOperations::getInstance()->print2d(R, M, N);
    MatrixOperations::getInstance()->print2d(Q, M, M);
    MatrixOperations::getInstance()->dmul(Q, R, X, M, N, M);
    MatrixOperations::getInstance()->print2d(X, M, N);


    return 0;
}