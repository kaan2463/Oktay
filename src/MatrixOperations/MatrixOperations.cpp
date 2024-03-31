#include "Exception.h"
#include "MatrixOperations.h"
#include <iostream>

#define DE 1.0e-8

MatrixOperations* MatrixOperations::INSTANCE = 0;

MatrixOperations* MatrixOperations::getInstance()
{
    if(INSTANCE == 0)
    {
        INSTANCE = new MatrixOperations;
    }
    return INSTANCE;
}

void MatrixOperations::axpy(double alpha, double* A, double beta, double* C, size_t M)
{
    for(size_t i = 0; i < M; i++)
    {
        C[i] = alpha * A[i] + beta * C[i];
    }
}

void MatrixOperations::dhad(double* A, double* B, double* C, size_t M)
{
    for(size_t i = 0; i < M; i++)
    {
        C[i] = A[i] * B[i];
    }
}

void MatrixOperations::dmul(double* A, double* B, double* C, size_t M, size_t N, size_t K)
{
    double value;
    size_t mk;
    size_t kn;
    size_t mn;
    for(size_t mi = 0; mi < M; mi++)
    {
        for(size_t ni = 0; ni < N; ni++)
        {
            mn = mi * N + ni;
            value = 0;
            for(size_t ki = 0; ki < K; ki++)
            {
                kn = ki * M + ni;
                mk = mi * K + ki;

                value += A[mk] * B[kn];
            }
            C[mn] = value;
        }
    }
}

void MatrixOperations::dsub(double* A, double* B, double* C, size_t M, size_t N)
{
    size_t mn;
    for(size_t mi = 0; mi < M; mi++)
    {
        for(size_t ni = 0; ni < N; ni++)
        {
            mn = mi * N + ni;

            C[mn] = A[mn] - B[mn];
        }
    }
}

void MatrixOperations::dadd(double* A, double* B, double* C, size_t M, size_t N)
{
    size_t mn;
    for(size_t mi = 0; mi < M; mi++)
    {
        for(size_t ni = 0; ni < N; ni++)
        {
            mn = mi * N + ni;

            C[mn] = A[mn] + B[mn];
        }
    }
}

void MatrixOperations::dcpy(double* A, double* B, size_t M, size_t N)
{
    size_t mn;
    for(size_t mi = 0; mi < M; mi++)
    {
        for(size_t ni = 0; ni < N; ni++)
        {
            mn = mi * N + ni;

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

void MatrixOperations::transpose2d(double* A, size_t M, size_t N)
{
    while(M > 1 && N > 1)
    {
        for(size_t i = 1; i < M; i++)
        {
            std::rotate(A + i, A + i * N, A + i * N + 1);
        }

        A += M;
        N -= 1;
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

double MatrixOperations::det2d(double* A, size_t M)
{
    double* L = new double[M * M];
    double* U = new double[M * M];

    LUDecomposition2d(A, L, U, M);

    double det = 1.0;
    size_t ii;
    for(size_t mi = 0; mi < M; mi++)
    {
        ii = mi * M + mi;
        det *= L[ii] * U[ii];
    }

    delete[] L;
    delete[] U;

    return det;
}

//inline double abs(double x)
//{
//    // breaks strict aliasing, but compiler writer knows this behavior for the platform
//    uint64_t i = reinterpret_cast<const std::uint64_t&>(x);
//    i &= 0x7FFFFFFFFFFFFFFFULL; // clear sign bit
//
//    return reinterpret_cast<const double&>(i);
//}

void MatrixOperations::inverse(double* A, double* B, size_t M)
{
    memset(B, 0, sizeof(double) * M * M);

    size_t ii, ij, jj, jk, ji, ik;
    double ci, cj;

    double* L = new double[M * M];

    memcpy(L, A, sizeof(double) * M * M);

    for(size_t i = 0; i < M; i++)
    {
        ii = i * M + i;
        B[ii] = 1.0;
    }

    for(size_t i = 0; i < M; i++)
    {
        MatrixOperations::print2d(L, M, M);
        MatrixOperations::print2d(B, M, M);
        ii = i * M + i;

        if(abs(L[ii]) < DE) // or equal 0.0 
        {
            for(size_t j = i + 1; j < M; j++)
            {
                ji = j * M + i;
                if(L[ji] > DE)
                {
                    for(size_t k = 0; k < M; k++)
                    {
                        ik = i * M + k;
                        jk = j * M + k;
                        swap(&L[ik], &L[jk]);
                        swap(&B[ik], &B[jk]);
                    }
                    break;
                }
            }
        }

        if(abs(L[ii]) < DE)
        {
            THROW_EXCEPTION("Singular Matrix\n");
            return;
        }

        ci = 1.0 / L[ii];
        for(size_t j = 0; j < M; j++)
        {
            ij = i * M + j;
            L[ij] = L[ij] * ci;
            B[ij] = B[ij] * ci;
        }
        MatrixOperations::print2d(L, M, M);
        MatrixOperations::print2d(B, M, M);

        for(size_t j = i + 1; j < M; j++)
        {
            ji = j * M + i;

            cj = L[ji];
            for(size_t k = 0; k < M; k++)
            {
                ik = i * M + k;
                jk = j * M + k;
                L[jk] -= L[ik] * cj;
                B[jk] -= B[ik] * cj;
            }
        }
    }

    for(size_t i = M - 1; i > 0; i--)
    {

        for(size_t j = 0; j < i; j++)
        {
            ji = (i - 1 - j) * M + i;
            cj = L[ji];
            for(size_t k = 0; k < M; k++)
            {
                ik = i * M + k;
                jk = (i - 1 - j) * M + k;
                L[jk] -= cj * L[ik];
                B[jk] -= cj * B[ik];
            }
        }
    }

}


void MatrixOperations::print2d(double* A, size_t M, size_t N)
{
    for(size_t mi = 0; mi < M; mi++)
    {
        for(size_t ni = 0; ni < N; ni++)
        {
            printf("%lf ", A[mi * N + ni]);
        }
        printf("\n");
    }
    printf("\n");
}