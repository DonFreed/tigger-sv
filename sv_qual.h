#ifndef SVQUAL
#define SVQUAL

#include <stdint.h>
#include "htslib/sam.h"
#include "plp2sv.h"

typedef struct {
    uint32_t qlen; // read (unitig) length
    float div; // divergence between the read and the reference allele
    uint32_t dp:8, qual:8, allele:8, is_fwd:1, tmp:7; // depth (quality), mapping quality, supported allele (0 for ref), is the read on the fwd strand
} read_qual_t;

typedef struct {
     int n, max;
     read_qual_t *a;
} qual_vec_t;

int get_qual_data(bam_hdr_t *h, int tid, int pos, int n, int *n_plp,const bam_pileup1_t **plp, khash_t(sv_hash) *sv_h, qual_vec_t *quals); 

#endif
