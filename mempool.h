#ifndef MEMPOLT
#define MEMPOLT

#include <stdint.h>
#include <stdlib.h>

/*******************
 *** Memory Pool ***
 *******************/

/* Allocate only memory pool for holding alleles and read support data */
/*  Substantially similar to the memory pool from ropebwt2/rope.c
    The license for ropebwt is included here
*/

/*
The MIT License

Copyright (c) 2013-2015 Broad Institute

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*/

#define MP_CHUNK_SIZE 0x10000 // 64 MB per chunk

typedef struct {
    int i, n_elems;
    int64_t top, max;
    uint8_t **mem;
} mempool_t;

mempool_t *mp_init();
void mp_destroy(mempool_t *mp);
inline void *mp_alloc(mempool_t *mp, size_t size);

#endif
