#ifndef OKTAY_MATRIX_OPERATIONS
#define OKTAY_MATRIX_OPERATIONS

#define MAX_ITER_EIG        255
#define EIG_DE              1.0e-6

class MatrixOperations
{
private:
    static MatrixOperations* INSTANCE;
    MatrixOperations() {}
public:
    static MatrixOperations* getInstance();
    ~MatrixOperations() {}

    /*
    * scalar multiply and sum
    * C = alpha * A + beta * C
    */
    void axpy(double alpha, double* A, double beta, double* C, size_t M);

    /*
    * Hadamard product
    * C_{i} = A_{i} * B_{i}
    */
    void dhad(double* A, double* B, double* C, size_t M);

    /*
    * Dot product
    * result = sum_{m} A_{m} * B_{m}
    */
    double ddot(double* A, double* B, size_t M);

    /*
    * Matrix Multiplication
    * C_{mn} = sum_{K} A_{mk} * B_{kn}
    */
    void dmul(double* A, double* B, double* C, size_t M, size_t N, size_t K);

    /*
    * Matrix Multiplication with alpha
    * C_{mn} = sum_{K} alpha * A_{mk} * B_{kn}
    */
    void dmul(double alpha, double* A, double* B, double* C, size_t M, size_t N, size_t K);

    /*
    * Matrix subtraction
    * C_{ij} = A_{ij} - B_{ij}
    */
    void dsub(double* A, double* B, double* C, size_t sz);

    /*
    * Matrix addition
    * C_{ij} = A_{ij} + B_{ij}
    */
    void dadd(double* A, double* B, double* C, size_t M, size_t N);

    /*
    * Matrix copy
    */
    void dcpy(double* dest, double* src, size_t sz);

    /*
    * In-Place transpose
    * A_{ij} = A_{ji}
    */
    void transpose2d(double* A, size_t M, size_t N);

    /*
    * LU decomposition of matrix A
    * Doolittle algorithm
    * L : letf triangular matrix
    * U : upper triangular matrix
    * A = LU
    */
    void lu(double* A, double* L, double* U, size_t M);

    /*
    * QR decomposition (Householder)
    * A = QR
    * Where Q is an orthogonal matrix (Q^{-1} = Q^{T})
    * And R is an upper triangular matrix
    * M >= N
    */
    void qr(double* A, double* Q, double* R, size_t M, size_t N);

    /*
    * Eigenvalue decomposition
    * QR factorization algorithm
    * E_{n}   : eigenvlaues
    * V_{n*n} : eigenvectors
    * A*V = E*V
    * A = V*E*V'
    */
    void eig(double* A, double* E, double* V, size_t M);

    /*
    * Singular Value Decomposition
    * A = U*E*V
    */
    void svd(double* A, double* U, double* E, double* V, size_t M, size_t N);


    /*
    * Determinant of square matrix A
    * return det(A)
    */
    double det2d(double* A, size_t M);

    /*
    * Inverse of Matrix
    * Gaussian Elimination Algorithm
    * C_{ij} = A^{-1}_{ij}
    */
    void inverse(double* A, double* C, size_t M);

    /*
    * Helper function for print Matrix
    */
    void print1d(double* A, size_t M);

    /*
    * Helper function for print Matrix
    */
    void print2d(double* A, size_t M, size_t N);
};
#endif