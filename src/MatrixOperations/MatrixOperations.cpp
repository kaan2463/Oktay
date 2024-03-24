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

void MatrixOperations::dsub(double* A, double* B, double* C, size_t m, size_t n)
{
    size_t mn;
    for(size_t mi = 0; mi < m; mi++)
    {
        for(size_t ni = 0; ni < n; ni++)
        {
            mn = mi * n + ni;

            C[mn] = A[mn] - B[mn];
        }
    }
}

void MatrixOperations::dadd(double* A, double* B, double* C, size_t m, size_t n)
{
    size_t mn;
    for(size_t mi = 0; mi < m; mi++)
    {
        for(size_t ni = 0; ni < n; ni++)
        {
            mn = mi * n + ni;

            C[mn] = A[mn] + B[mn];
        }
    }
}

void MatrixOperations::dcpy(double* A, double* B, size_t m, size_t n)
{
    //memcpy(A, B, m * n * sizeof(double));

    size_t mn;
    for(size_t mi = 0; mi < m; mi++)
    {
        for(size_t ni = 0; ni < n; ni++)
        {
            mn = mi * n + ni;

            A[mn] = B[mn];
        }
    }
}


//Helper swap
inline void swap(double* a, double* b)
{
    double xa = *a;
    *a = *b;
    *b = xa;

}

//Helper rotate
inline void rotate(double* first, double* middle, double* last)
{
    double* next = middle;
    while(first != next)
    {
        swap(first++, next++);
        if(next == last) next = middle;
        else if(first == middle) middle = next;
    }
}

void MatrixOperations::transpose2d(double* A, size_t m, size_t n)
{
    while(m > 1 && n > 1)
    {
        for(size_t i = 1; i < m; i++)
        {
            std::rotate(A + i, A + i * n, A + i * n + 1);
        }

        A += m;
        n -= 1;
    }
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