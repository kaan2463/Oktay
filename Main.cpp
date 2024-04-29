#include <MathematicalOperations.h>
#include <Matrix.h>

#include <iostream>
#include <MatrixOperations.h>
using namespace std;

#define MAT MatrixOperations::getInstance()


int main()
{

    size_t M = 4;

    //double* A = new double[M * N];
   // V = [9 - 4 - 2 0; -56 32 - 28 44; -14 - 14 6 - 14; 42 - 33 21 - 45]
    double A[] = { 9, -4, -2, 0,
        -56, 32 ,-28, 44,
        -14, -14, 6, -14,
        42, -33, 21, -45 };
    double B[] = {
        2,0,0,
        0,3,4,
        0,4,9 };


    double* E = new double[M * M];
    double* V = new double[M * M];
    double* IV = new double[M * M];

    MAT->eig(A, E, V, M);
    MAT->inverse(V, IV, M);
    MAT->print1d(E, M);
    MAT->print2d(V, M, M);
    MAT->print2d(IV, M, M);

    double* T = new double[M * M];
    double* C = new double[M * M];
    MAT->dmul(IV, A, C, M, M, M);
    MAT->dmul(C, V, T, M, M, M);
    MAT->print2d(T, M, M);
    MAT->eig(T, E, V, M);
    MAT->print1d(E, M);
    MAT->print2d(V, M, M);

    return 0;
}