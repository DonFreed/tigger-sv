#include <stdio.h>
#include "mempool.h"

mempool_t *mp_init()
{
    mempool_t *mp;
    mp = calloc(1, sizeof(mempool_t));
    mp->i = mp->n_elems = MP_CHUNK_SIZE;
    mp->top = -1;
    return mp;
}

void mp_destroy(mempool_t *mp)
{
    int64_t i;
    for (i = 0; i <= mp->top; ++i) { free(mp->mem[i]); }
    free(mp->mem);
    free(mp);
}

inline void *mp_alloc(mempool_t *mp, size_t size)
{
    int64_t i;
    if (mp->i + size >= mp->n_elems) {
        if (++mp->top == mp->max) {
            mp->max = mp->max ? mp->max << 1 : 1;
            mp->mem = realloc(mp->mem, sizeof(void*) * mp->max);
        }
        mp->mem[mp->top] = calloc(mp->n_elems, sizeof(uint8_t));
        mp->i = 0;
    }
    i = mp->i;
    mp->i += size;
    return mp->mem[mp->top] + i;
}
