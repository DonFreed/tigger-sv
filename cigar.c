#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include "htslib/sam.h"
#include "cigar.h"

inline void parse_cigar(uint8_t *qual, int32_t l_qseq, const uint32_t *cigar, int n_cigar, cigar_res_t *res)
{
    int k;
    if (qual) {
        res->clip_q[0] = qual[0]; 
        res->clip_q[1] = qual[l_qseq];
    }
    res->clip[0] = 0; res->clip[1] = 0; res->qlen = 0; res->rlen = 0; res->ins = 0; res->del = 0;
    for (k = 0; k < n_cigar; ++k) {
        int op = bam_cigar_op(cigar[k]);
        int oplen = bam_cigar_oplen(cigar[k]);
        if (oplen == 0) continue;
        if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP) res->qlen += oplen;
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            res->clip[!!k] = oplen;
            if (qual && op == BAM_CSOFT_CLIP) {
                res->clip_q[!!k] = !k? qual[oplen] : qual[res->clip[0] + res->qlen - 1];
            }
        }
        if (bam_cigar_type(op)&2) res->rlen += oplen;
        if (op == BAM_CINS) res->ins += oplen;
        else if (op == BAM_CDEL) res->del += oplen;
    }
    return;
}

inline void str2cigar(const char *s, i32array_t *cigar, int *n_cigar)
{
    const char *_s;
    uint32_t l = 0;
    uint8_t op = 0;
    cigar->n = 0;
    for (_s = s; *_s; ++_s) {
        if (isdigit(*_s)) {
            l *= 10;
            l += *_s - '0';
        } else {
            if (*_s == 'M') { op = 0; }
            else if (*_s == 'I') { op = 1; }
            else if (*_s == 'D') { op = 2; }
            else if (*_s == 'N') { op = 3; }
            else if (*_s == 'S') { op = 4; }
            else if (*_s == 'H') { op = 5; }
            else if (*_s == 'P') { op = 6; }
            else if (*_s == '=') { op = 7; }
            else if (*_s == 'X') { op = 8; }
            else if (*_s == 'B') { op = 9; }
            else { cigar->a = 0; *n_cigar = -1; return; }
            
            if (cigar->n >= cigar->max) {
                cigar->max = cigar->n ? cigar->n + 1 : 32;
                kroundup32(cigar->max);
                if (!(cigar->a = (uint32_t*)realloc(cigar->a, cigar->max * sizeof(uint32_t)))) {
                    cigar->a = 0; *n_cigar = -1; return;
                }
            }
            cigar->a[cigar->n++] = (l << 4) | op;
            l = 0;
        }
    }
    *n_cigar = cigar->n;
}

