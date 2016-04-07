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
} cigar_res_t;

inline void parse_cigar(const uint32_t *cigar, int n_cigar, cigar_res_t *res);

inline void str2cigar(const char *s, i32array_t *cigar, int *n_cigar);

inline int parse_sa_tag(bam_hdr_t *h, kstring_t *sa, int is_front, int is_rev, int sv_qbeg, cigar_res_t *res, int *sa_is_rev, int *sa_qbeg, int *sa_tid, int *sa_pos);

inline uint32_t cigar_get_qual(uint8_t *qual, int tdist, const uint32_t *cigar, int n_cigar);

#endif
