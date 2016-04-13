#ifndef TGGRGNO
#define TGGRGNO

#include "htslib/sam.h"
#include "sv_qual.h"

void print_header(bam_hdr_t *h, int optind, int n, char *argv[]);

int genotype_sv(bam_hdr_t *h, int n, khash_t(sv_geno) *geno_h);

#endif
