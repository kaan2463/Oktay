#ifndef OKTAY_MathematicalOperationsSTL
#define OKTAY_MathematicalOperationsSTL

#include <MathematicalOperations.h>

class MathematicalOperationsSTL : public MathematicalOperations
{
    static MathematicalOperationsSTL* INSTANCE;
    MathematicalOperationsSTL() {}
public:
    static MathematicalOperations* getInstance();
    ~MathematicalOperationsSTL() {}

    double abs(double x);

    double sqrt(double x);

    double pow10(double x);

    double powInv10(double x);

    double pow(double x, double n);

    double exp(double x);

    double sin(double x);

    double cos(double x);

    double tan(double x);

    double arcsin(double x);

    double arccos(double x);

    double arctan(double x);

    double ln(double x);

    double log(double base, double x);

    double log2(double x);

    double log10(double x);
};
#endif
