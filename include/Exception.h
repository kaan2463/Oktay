#ifndef OKTAY_EXCEPTION_H
#define OKTAY_EXCEPTION_H

#include <io.h>

class RuntimeException
{
public:
    RuntimeException(const char* msg)
    {
        printfW("ERROR : %s \n", msg);
    }
};

#endif

