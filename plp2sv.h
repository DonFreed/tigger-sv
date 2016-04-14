#ifndef PLP2SV
#define PLP2SV

#include <stdint.h>
#include "htslib/khash.h"
#include "htslib/sam.h"

typedef struct {
    uint64_t id; // position of bp1 << 32 | position of bp2 where 
                 //  pos bp1 is defined as min(pos_bp1, pos_bp2)
    int allele; // the numeric id of the allele
    int32_t tid1;
    int32_t pos1;
    int ori1; // 0 = read extends to the right of bp; 1 = read extends to the left of bp.
    int32_t tid2;
    int32_t pos2;
    int ori2;
    char type;
    int qdist; // size of the SV on the query (unitig)
} sv_t;

typedef struct {
    int n, max;
    sv_t *sv;
} sv_vec_t;

KHASH_MAP_INIT_INT64(sv_hash, sv_t)

/*
 * plp2sv
 *  tid is the integer id of the reference chromosome (input)
 *  pos is the current genomic position (input)
 *  n is the number of individuals (intput)
 *  n_plp is the number of reads in each individual (input)
 *  plp is a pileup of reads (input)
 *  sv is data on detected structural variants (output)
 *  return is the number of detected structural variants.
 */
int plp2sv(bam_hdr_t *h, int tid, int pos, int n, int *n_plp, const bam_pileup1_t **plp, khash_t(sv_hash) *sv_h);

#endif
