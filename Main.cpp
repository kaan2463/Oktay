#include <MatrixOperations.h>

#define MAT MatrixOperations::getInstance()



int main()
{
    const size_t M = 3;
    double C[M * M];
    double D[M * M];
    double A[] = {
        1,2,3,
        123,54,22,
        2,4,6
    };

    MAT->inverse(A, C, M); //throw exception because A is singular

    MAT->print2d(C, M, M);
    MAT->dmul(A, C, D, M, M, M);
    MAT->print2d(D, M, M);
}