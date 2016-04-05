#ifndef TGGRARR
#define TGGRARR

#include <stdint.h>

typedef struct {
    int n, max, *a;
} iarray_t;

typedef struct {
    int n, max;
    uint32_t *a;
} i32array_t;

#endif
