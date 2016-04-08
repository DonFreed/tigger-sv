#ifndef SUMQUAL
#define SUMQUAL

#include <stdint.h>
#include "htslib/khash.h"
#include "mempool.h"
#include "plp2sv.h"
#include "sv_qual.h"

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t genotyped:1, ori1:1, ori2:1, tmp:29; // whether the allele has been genotyped, orientation of the SV at the site. Same as ori in plp2sv.h
} allele_t; 

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t mq_sum; // sums of squares of mapq
    double div_sum;
    uint32_t qlen_sum;
    uint32_t n_reads; // number of reads found at the site
    uint8_t n_alleles; // number of identified alleles
    uint16_t *read_data; // pointer to an array of (n * a) with n = # of individuals and a = # of alleles
    allele_t *alleles; // pointer to an array of allele_t objects for each detected variant allele
} qual_sum_t;

KHASH_MAP_INIT_INT64(sv_geno, qual_sum_t);

int summarize_qual(uint32_t tid, uint32_t pos, int n, int n_alleles, khash_t(sv_hash) *sv_h, qual_vec_t *quals, mempool_t *mp, khash_t(sv_geno) *geno_h);

#endif
