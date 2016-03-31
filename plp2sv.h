#ifndef PLP2SV
#define PLP2SV

#include <stdint.h>

typedef struct {
    uint64_t id; // position of bp1 << 32 | position of bp2 where 
                 //  pos bp1 is defined as min(pos_bp1, pos_bp2)
    int32_t tid;
    int32_t pos;
    int32_t to_pos; // position of the other alignment of the breakpoint
    int qlen; // length of the alignment along the read
    int rlen; // length of the alignment along the reference
    int flag; // sam flag
    int mapq;
    int qbeg; // how far along the read until the alignment begins
    int clip[2]; // size of soft or hard clips
    int ins; // the length of insertions in the alignment
    int del; // the length of deletions in the alignment
    int nm; // the number of mismatches in the alignment
    int score; // the alignment score
    int tipq; // the depth at the breakpoint
    double diff; // the rate of mismatches along the reference
} sv_t;

typedef struct {
    int n, max;
    sv_t *sv;
} sv_vec_t;

typedef struct {
    int n, max, *a;
} iarray_t;

typedef struct {
    int n, max;
    uint32_t *a;
} i32array_t;

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
inline int plp2sv(int tid, int pos, int n; int *n_plp; bam_pileup1_t *plp; sv_vec_t *sv)

#endif
