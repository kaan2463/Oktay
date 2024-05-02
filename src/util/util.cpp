#include <util.h>

#include <string.h>

void* memsetW(void* ptr, int value, size_t num)
{
    return memset(ptr, value, num);
}

void* memcpyW(void* dest, const void* src, size_t count)
{
    return memcpy(dest, src, count);
}

