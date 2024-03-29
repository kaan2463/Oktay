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
    void axpy(double alpha, double* A, double beta, double* C, size_t m);

    /*
    * Hadamard product
    * C_{i} = A_{i} * B_{i}
    */
    void dhad(double* A, double* B, double* C, size_t m);

    /*
    * Matrix Multiplication
    * C_{mn} = sum_{k} A_{mk} * B_{kn}
    */
    void dmul(double* A, double* B, double* C, size_t m, size_t n, size_t k);

    /*
    * Matrix subtraction
    * C_{ij} = A_{ij} - B_{ij}
    */
    void dsub(double* A, double* B, double* C, size_t m, size_t n);

    /*
    * Matrix addition
    * C_{ij} = A_{ij} + B_{ij}
    */
    void dadd(double* A, double* B, double* C, size_t m, size_t n);

    /*
    * Matrix copy
    * B_{ij} = A_{ij}
    */
    void dcpy(double* A, double* B, size_t m, size_t n);

    /*
    * In-Place transpose
    * A_{ij} = A_{ji}
    */
    void transpose2d(double* A, size_t m, size_t n);

    /*
    * Helper function for print Matrix
    */
    void print2d(double* A, size_t m, size_t n);
};
#endif