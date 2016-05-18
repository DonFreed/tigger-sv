#ifndef TGGSEG
#define TGGSEG

#include "htslib/khash.h"
#include "genotype.h"
#include "ped.h"

double log_segregation(genotype_t *gts, int n, double mi_prob, khash_t(ped) *ped_h);

#endif
