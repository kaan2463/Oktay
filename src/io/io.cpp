#include <stdarg.h>
#include <stdio.h>

#include <io.h>

void printfW(const char* format, ...)
{
    va_list args;
    va_start(args, format);

    vprintf(format, args);

    va_end(args);
}