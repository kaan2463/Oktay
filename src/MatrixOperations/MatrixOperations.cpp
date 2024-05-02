#include <Exception.h>
#include <io.h>
#include <MathematicalOperations.h>
#include <MatrixOperations.h>
#include <util.h>

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

double MatrixOperations::ddot(double* A, double* B, size_t M)
{
    double sum = 0.0;
    for(size_t i = 0; i < M; i++)
    {
        sum += A[i] * B[i];
    }
    return sum;
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
                kn = ki * N + ni;
                mk = mi * K + ki;

                value += A[mk] * B[kn];
            }
            C[mn] = value;
        }
    }
}

void MatrixOperations::dmul(double alpha, double* A, double* B, double* C, size_t M, size_t N, size_t K)
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
                kn = ki * N + ni;
                mk = mi * K + ki;

                value += A[mk] * B[kn];
            }
            C[mn] = alpha * value;
        }
    }
}

void MatrixOperations::dsub(double* A, double* B, double* C, size_t sz)
{
    for(size_t i = 0; i < sz; i++)
    {
        C[i] = A[i] - B[i];
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

void MatrixOperations::dcpy(double* dest, double* src, size_t sz)
{
    for(size_t i = 0; i < sz; i++)
    {
        dest[i] = src[i];
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
            rotate(A + i, A + i * N, A + i * N + 1);
        }

        A += M;
        N -= 1;
    }
}

void MatrixOperations::lu(double* A, double* L, double* U, size_t M)
{
    size_t ik, ii, kj, km, mj, kk, im, mk;

    memsetW(L, 0, sizeof(double) * M * M);
    memsetW(U, 0, sizeof(double) * M * M);

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

/*
* Generate M*M Identity Matrix
*/
inline static void eye(double* A, size_t M)
{
    size_t mn;
    for(size_t m = 0; m < M; m++)
    {
        for(size_t n = 0; n < M; n++)
        {
            mn = m * M + n;
            if(m == n)
            {
                A[mn] = 1.0;
            }
            else
            {
                A[mn] = 0.0;
            }

        }
    }
}

inline static double norm(double* V, size_t M)
{
    double sum = 0.0;
    for(size_t i = 0; i < M; i++)
    {
        sum += V[i];
    }
    return MathematicalOperations::getInstance()->sqrt(sum);
}

inline static double norm(double* V, size_t stride, size_t M)
{
    double sum = 0.0;
    for(size_t i = 0; i < M; i++)
    {
        sum += V[i * stride] * V[i * stride];
    }
    return MathematicalOperations::getInstance()->sqrt(sum);
}

//implementation: https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
void MatrixOperations::qr(double* A, double* Q, double* R, size_t M, size_t N)
{
    if(M < N)
    {
        THROW_EXCEPTION("M can not be less than N!");
        return;
    }
    eye(Q, M);
    dcpy(R, A, M * N);

    double normx, s, u1, tau;
    double t;
    double* C = new double[M * M];
    double* T = new double[M * M];
    double* w = new double[M];
    size_t ii;

    for(size_t i = 0; i < N; i++)
    {
        ii = i * N + i;
        normx = norm(&R[i + i * N], N, M - i);
        s = R[ii] < 0 ? 1.0 : -1.0;
        u1 = R[ii] - s * normx;
        w[0] = 1.0;
        for(size_t j = 1; j < M - i; j++)
        {
            w[j] = R[i + (j + i) * N] / u1;
        }
        tau = -s * u1 / normx;

        dmul(tau, w, w, C, M - i, M - i, 1); // C = w*w' => C=C'
        dmul(C, &R[i * N], T, M - i, N, M - i);
        dsub(&R[i * N], T, &R[i * N], (M - i) * N);

        for(size_t j = 0; j < M; j++)
        {
            for(size_t k = 0; k < M - i; k++)
            {
                t = 0;
                for(size_t l = 0; l < M - i; l++)
                {
                    t += Q[i + j * M + l] * C[l * (M - i) + k];
                }
                T[j * (M - i) + k] = t;
            }
        }

        for(size_t j = 0; j < M; j++)
        {
            for(size_t k = 0; k < M - i; k++)
            {
                Q[i + j * M + k] -= T[j * (M - i) + k];
            }
        }

    }
    delete[] C;
    delete[] T;
    delete[] w;

}

void MatrixOperations::eig(double* A, double* E, double* V, size_t M)
{
    double* A0 = new double[M * M];
    double* Q = new double[M * M];
    double* R = new double[M * M];
    double* T = new double[M * M];



    size_t iter = 0;
    size_t indexL;
    double sumL;
    dcpy(A0, A, M * M);
    eye(V, M);
    do
    {
        qr(A0, Q, R, M, M);
        dmul(R, Q, A0, M, M, M);
        dmul(Q, V, T, M, M, M);
        dcpy(V, T, M * M);

        sumL = 0.0;
        for(size_t i = 0; i < M - 1; i++)
        {
            for(size_t j = 0; j < i + 1; j++)
            {
                indexL = i * M + j + M;
                sumL += MathematicalOperations::getInstance()->abs(A0[indexL]);
            }
        }

        iter++;
    } while(EIG_DE < MathematicalOperations::getInstance()->abs(sumL) && iter < MAX_ITER_EIG);

    for(size_t i = 0; i < M; i++)
    {
        E[i] = A0[i * M + i];
    }

    delete[] A0;
    delete[] Q;
    delete[] R;
    delete[] T;
}

void MatrixOperations::svd(double* A, double* U, double* E, double* V, size_t M, size_t N)
{
    size_t P = M > N ? M : N;
    double* T = new double[P * P];
    double* AT = new double[M * N];
    double* EE = new double[P];

    dcpy(AT, A, M * N);
    transpose2d(AT, M, N);

    dmul(A, AT, T, M, M, N);

    eig(T, EE, U, M);

    if(M == P)
    {
        for(size_t i = 0; i < M; i++)
        {
            E[i] = MathematicalOperations::getInstance()->sqrt(EE[i]);
        }
    }

    dmul(AT, A, T, N, N, M);

    eig(T, EE, V, N);

    if(N == P)
    {
        for(size_t i = 0; i < N; i++)
        {
            E[i] = MathematicalOperations::getInstance()->sqrt(EE[i]);
        }
    }


    delete[] T;
    delete[] AT;
    delete[] EE;

}

double MatrixOperations::det2d(double* A, size_t M)
{
    double* L = new double[M * M];
    double* U = new double[M * M];

    lu(A, L, U, M);

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
    memsetW(B, 0, sizeof(double) * M * M);

    size_t ii, ij, jk, ji, ik;
    double ci, cj;

    double* L = new double[M * M];

    memcpyW(L, A, sizeof(double) * M * M);

    for(size_t i = 0; i < M; i++)
    {
        ii = i * M + i;
        B[ii] = 1.0;
    }

    for(size_t i = 0; i < M; i++)
    {
        ii = i * M + i;

        if(MathematicalOperations::getInstance()->abs(L[ii]) < DE) // or equal 0.0 
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

        if(MathematicalOperations::getInstance()->abs(L[ii]) < DE)
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

void MatrixOperations::print1d(double* A, size_t M)
{
    for(size_t mi = 0; mi < M; mi++)
    {
        printfW("%lf ", A[mi]);
    }
    printfW("\n");
    printfW("\n");
}


void MatrixOperations::print2d(double* A, size_t M, size_t N)
{
    for(size_t mi = 0; mi < M; mi++)
    {
        for(size_t ni = 0; ni < N; ni++)
        {
            printfW("%lf ", A[mi * N + ni]);
        }
        printfW("\n");
    }
    printfW("\n");
}