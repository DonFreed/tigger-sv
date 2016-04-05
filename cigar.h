#ifndef TIGCIG
#define TIGCIG

#include <stdint.h>
#include "array.h"

typedef struct {
    int clip[2]; // size of soft or hard clips
    int qlen; // length of the alignment along the read
    int rlen; // length of the alignment along the reference
    int ins; // the length of insertions in the alignment
    int del; // the length of deletions in the alignment
    int clip_q[2]; // the quality score (depth) at the tips of the alignment
} cigar_res_t;

inline void parse_cigar(uint8_t *qual, int32_t l_qseq, const uint32_t *cigar, int n_cigar, cigar_res_t *res);

inline void str2cigar(const char *s, i32array_t *cigar, int *n_cigar);

#endif
