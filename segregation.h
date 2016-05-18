#ifndef TGGSEG
#define TGGSEG

#include "htslib/khash.h"
#include "genotype.h"
#include "ped.h"

void gt_to_alleles(int gt, int *alleles);
double log_segregation(genotype_t *gts, int n, double mi_prob, khash_t(ped) *ped_h);
int n_mi(genotype_t *gts, int n, khash_t(ped) *ped_h);

#endif
