#ifndef TIGPED
#define TIGPED

#include "htslib/khash.h"

typedef struct {
    int fa_col;
    int mo_col;
} parents_t;

KHASH_MAP_INIT_STR(colmap, int);
KHASH_MAP_INIT_INT(ped, parents_t);

khash_t(colmap) *map_samples(char **samples, int n);
khash_t(ped) *read_ped(const char *fnped, khash_t(colmap) *smp_cols);

#endif
