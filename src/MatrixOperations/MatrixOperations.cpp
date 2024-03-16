#include "MatrixOperations.h"
#include <iostream>

MatrixOperations* MatrixOperations::INSTANCE = 0;

MatrixOperations* MatrixOperations::getInstance()
{
    if(INSTANCE == 0)
    {
        INSTANCE = new MatrixOperations;
    }
    return INSTANCE;
}



void MatrixOperations::axpy(double alpha, double* A, double beta, double* C, size_t m)
{
    for(size_t i = 0; i < m; i++)
    {
        C[i] = alpha * A[i] + beta * C[i];
    }
}

void MatrixOperations::dhad(double* A, double* B, double* C, size_t m)
{
    for(size_t i = 0; i < m; i++)
    {
        C[i] = A[i] * B[i];
    }
}

void MatrixOperations::dmul(double* A, double* B, double* C, size_t m, size_t n, size_t k)
{
    double value;
    size_t mk = 0;
    size_t km = 0;
    size_t mn = 0;
    for(size_t mi = 0; mi < m; mi++)
    {
        for(size_t ni = 0; ni < n; ni++)
        {
            mn = mi * n + ni;
            value = 0;
            for(size_t ki = 0; ki < k; ki++)
            {
                km = ki * m + mi;
                mk = mi * k + ki;

                value += A[mk] * B[km];
            }
            C[mn] = value;
        }
    }
}

// TODO kaan: Implement In-Place transpose algorithm! 
void MatrixOperations::transpose2d(double* A, size_t m, size_t n)
{

    // size_t mn = 0;
    // size_t nm = 0;
    // double value;
    // 
    // for(size_t mi = 0; mi < m; mi++)
    // {
    //     for(size_t ni = 0; ni < mi + 1; ni++)
    //     {
    //         mn = mi * n + ni;
    //         nm = ni * m + mi;
    //         value = A[mn];
    //         A[mn] = A[nm];
    //         A[nm] = value;
    //     }
    // }

}


void MatrixOperations::print2d(double* A, size_t m, size_t n)
{
    for(size_t mi = 0; mi < m; mi++)
    {
        for(size_t ni = 0; ni < n; ni++)
        {
            printf("%lf ", A[mi * n + ni]);
        }
        printf("\n");
    }
    printf("\n");
}