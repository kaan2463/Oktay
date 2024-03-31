#ifndef OKTAY_MATRIX_OPERATIONS
#define OKTAY_MATRIX_OPERATIONS
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
    * Matrix Multiplication
    * C_{mn} = sum_{K} A_{mk} * B_{kn}
    */
    void dmul(double* A, double* B, double* C, size_t M, size_t N, size_t K);

    /*
    * Matrix subtraction
    * C_{ij} = A_{ij} - B_{ij}
    */
    void dsub(double* A, double* B, double* C, size_t M, size_t N);

    /*
    * Matrix addition
    * C_{ij} = A_{ij} + B_{ij}
    */
    void dadd(double* A, double* B, double* C, size_t M, size_t N);

    /*
    * Matrix copy
    * B_{ij} = A_{ij}
    */
    void dcpy(double* A, double* B, size_t M, size_t N);

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
    void LUDecomposition2d(double* A, double* L, double* U, size_t M);

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
    void print2d(double* A, size_t M, size_t N);
};
#endif