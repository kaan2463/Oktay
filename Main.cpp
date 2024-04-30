#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
#include <MatrixOperations.h>
using namespace std;

#define MAT MatrixOperations::getInstance()


int main()
{

    size_t M = 2;
    size_t N = 3;

    //double* A = new double[M * N];
   // V = [9 - 4 - 2 0; -56 32 - 28 44; -14 - 14 6 - 14; 42 - 33 21 - 45]
    double A[] = { 9, -4, -2, 0,
        -56, 32 ,-28, 44,
        -14, -14, 6, -14,
        42, -33, 21, -45 };
    double B[] = {
        2,0,
        0,3,
        0,4 };

    double C[] = { 3, 2, 2,
2, 3, -2 };




    double* E = new double[M > N ? M : N];
    double* V = new double[N * N];
    double* U = new double[M * M];


    MAT->svd(C, U, E, V, M, N);

    MAT->print2d(U, M, M);
    MAT->print1d(E, M > N ? M : N);
    MAT->print2d(V, N, N);

    return 0;
}