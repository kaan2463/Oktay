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
    size_t mk;
    size_t kn;
    size_t mn;
    for(size_t mi = 0; mi < m; mi++)
    {
        for(size_t ni = 0; ni < n; ni++)
        {
            mn = mi * n + ni;
            value = 0;
            for(size_t ki = 0; ki < k; ki++)
            {
                kn = ki * m + ni;
                mk = mi * k + ki;

                value += A[mk] * B[kn];
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

void MatrixOperations::LUDecomposition2d(double* A, double* L, double* U, size_t M)
{
    size_t ik, ii, kj, km, mj, kk, im, mk;

    memset(L, 0, sizeof(double) * M * M);
    memset(U, 0, sizeof(double) * M * M);

    for(size_t i = 0; i < M; ++i)
    {
        ii = i * M + i;
        U[i] = A[i];
        L[ii] = 1.0;
        L[i * M] = A[i * M] / U[0];
    }
    for(size_t k = 1; k < M; ++k)
    {
        kk = k * M + k;
        for(size_t j = k; j < M; ++j)
        {
            kj = k * M + j;
            U[kj] = A[kj];
            for(size_t m = 0; m < k; ++m)
            {
                km = k * M + m;
                mj = m * M + j;
                U[kj] -= L[km] * U[mj];
            }
        }
        for(size_t i = k + 1; i < M; ++i)
        {
            ik = i * M + k;
            L[ik] = A[ik];
            for(size_t m = 0; m < k; ++m)
            {
                im = i * M + m;
                mk = m * M + k;
                L[ik] -= L[im] * U[mk];
            }
            L[ik] /= U[kk];
        }
    }
}

double MatrixOperations::det2d(double* A, size_t m)
{
    double* L = new double[m * m];
    double* U = new double[m * m];

    LUDecomposition2d(A, L, U, m);

    double det = 1.0;
    size_t ii;
    for(size_t mi = 0; mi < m; mi++)
    {
        ii = mi * m + mi;
        det *= L[ii] * U[ii];
    }

    delete[] L;
    delete[] U;

    return det;
}

void MatrixOperations::inverse(double* A, double* B, size_t m)
{
    // TODO kaan: implement
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