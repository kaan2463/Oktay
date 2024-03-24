#ifndef OKTAY_MATRIX_H
#define OKTAY_MATRIX_H
class Matrix
{
private:
    size_t sz;
    double* data;
public:

    Matrix()
    {
        sz = 0;
        data = 0;
    }

    Matrix(double* data, size_t sz)
    {
        this->sz = sz;
        this->data = data;
    }

    Matrix(size_t sz)
    {
        this->sz = sz;
        data = new double[sz];
    }

    ~Matrix()
    {
        delete[] data;
    }

    size_t Size()
    {
        return sz;
    }

    void Size(size_t sz)
    {
        this->sz = sz;
    }

    double* Data()
    {
        return data;
    }

    void Data(double* data, size_t)
    {
        this->data = data;
    }

};

class Matrix2d : public Matrix
{
private:
    size_t m; // number of rows
    size_t n; // number of columns
public:

    Matrix2d(const Matrix2d& A);


    Matrix2d(size_t m, size_t n) : Matrix(m* n)
    {
        this->m = m;
        this->n = n;
    }

    Matrix2d(double* data, size_t m, size_t n) : Matrix(data, m* n)
    {
        this->m = m;
        this->n = n;
    }
    ~Matrix2d()
    {
    }

    void Dim2d(size_t m, size_t n)
    {
        this->m = m;
        this->n = n;
    }

    void print();

    // operators
    Matrix2d operator*(Matrix2d A);
    Matrix2d operator+(Matrix2d A);
    Matrix2d operator-(Matrix2d A);
    Matrix2d& operator=(Matrix2d A);
    Matrix2d& operator=(Matrix2d& A);
    Matrix2d operator+();
};
#endif