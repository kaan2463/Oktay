#include "Matrix.h"
#include "MatrixOperations.h"

Matrix2d::Matrix2d(const Matrix2d& A) :Matrix(m* n)
{
    this->m = A.m;
    this->n = A.n;
    MatrixOperations::getInstance()->dcpy(Data(), ((Matrix)A).Data(), m * n);
}

void Matrix2d::print()
{
    MatrixOperations::getInstance()->print2d(Data(), m, n);
}

Matrix2d Matrix2d::operator*(Matrix2d& A)
{
    Matrix2d C(this->m, A.n);
    MatrixOperations::getInstance()->dmul(Data(), A.Data(), C.Data(), this->m, A.n, A.m);
    return C;
}

Matrix2d Matrix2d::operator+(Matrix2d& A)
{
    Matrix2d C(A.m, A.n);
    MatrixOperations::getInstance()->dadd(Data(), A.Data(), C.Data(), A.m, A.n);
    return C;
}

Matrix2d Matrix2d::operator-(Matrix2d& A)
{
    Matrix2d C(A.m, A.n);
    MatrixOperations::getInstance()->dsub(Data(), A.Data(), C.Data(), A.m * A.n);
    return C;
}

Matrix2d& Matrix2d::operator=(Matrix2d& A)
{
    m = A.m;
    n = A.n;
    MatrixOperations::getInstance()->dcpy(Data(), A.Data(), A.m * A.n);
    return *this;
}

Matrix2d& Matrix2d::operator=(const Matrix2d& A)
{
    m = A.m;
    n = A.n;
    MatrixOperations::getInstance()->dcpy(Data(), ((Matrix)A).Data(), A.m * A.n);
    return *this;
}

Matrix2d Matrix2d::operator+()
{
    Matrix2d C(m, n);
    MatrixOperations::getInstance()->dcpy(C.Data(), Data(), m * n);
    MatrixOperations::getInstance()->transpose2d(C.Data(), m, n);
    C.Dim2d(n, m);
    return C;
}
