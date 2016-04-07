#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "cigar.h"

inline uint32_t cigar_get_qual(uint8_t *qual, int tdist, const uint32_t *cigar, int n_cigar)
{
    int k;
    int rdist = 0; // distance along the reference
    int qdist = 0; // distance along the query (read)
    for (k = 0; k < n_cigar; ++k) {
        int op = bam_cigar_op(cigar[k]);
        int oplen = bam_cigar_oplen(cigar[k]);
        if (oplen == 0) continue;
        if (bam_cigar_type(op)&2) {
            if (rdist + oplen > tdist) {
                if (bam_cigar_type(op)&1) {
                    return qual[qdist + (tdist - rdist)];
                } else {
                    return qual[qdist];
                }
            } else {
                rdist += oplen;
            }
        }
        if (bam_cigar_type(op)&1) {
            qdist += oplen;
        }
    }
    return UINT32_MAX;
}

inline void parse_cigar(const uint32_t *cigar, int n_cigar, cigar_res_t *res)
{
    int k;
    res->clip[0] = 0; res->clip[1] = 0; res->qlen = 0; res->rlen = 0; res->ins = 0; res->del = 0;
    for (k = 0; k < n_cigar; ++k) {
        int op = bam_cigar_op(cigar[k]);
        int oplen = bam_cigar_oplen(cigar[k]);
        if (oplen == 0) continue;
        if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP) res->qlen += oplen;
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            res->clip[!!k] = oplen;
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

inline int parse_sa_tag(bam_hdr_t *h, kstring_t *sa, int is_front, int is_rev, int sv_qbeg, cigar_res_t *res, int *sa_is_rev, int *sa_qbeg, int *sa_tid, int *sa_pos)
{
    iarray_t alns = {0, 0, 0};
    iarray_t fields = {0, 0, 0};
    i32array_t sa_cigar = {0, 0, 0};
    int k, sa_n_cigar, is_sv, sa_idx = -1, _sa_is_rev, _sa_qbeg;
    uint32_t sa_off[8] = {0};
    cigar_res_t _res;
    alns.n = ksplit_core(sa->s, ';', &alns.max, &alns.a);
    for (k = 0; k < alns.n; ++k) {
        fields.n = ksplit_core(sa->s + alns.a[k], ',', &fields.max, &fields.a);
        sa_cigar.n = 0;
        str2cigar(sa->s + alns.a[k] + fields.a[3], &sa_cigar, &sa_n_cigar);
        parse_cigar(sa_cigar.a, sa_n_cigar, &_res);
        _sa_is_rev = *(sa->s + alns.a[k] + fields.a[2]) == '-';
        _sa_qbeg = _sa_is_rev ? _res.clip[1] : _res.clip[0];
        is_sv = 0;
        if (is_front ^ (!is_rev)) {
            if (sv_qbeg < _sa_qbeg) {
                if (sa_idx >= 0) {
                    if (_sa_qbeg < *sa_qbeg) {
                        is_sv = 1;
                    }
                } else {
                    is_sv = 1;
                }
            }
        } else {
            if (sv_qbeg > _sa_qbeg) {
                if (sa_idx >= 0) {
                    if (_sa_qbeg > *sa_qbeg) {
                        is_sv = 1;
                    }
                } else {
                    is_sv = 1;
                }
            }
        }
        if (is_sv) {
            sa_idx = k;
            *sa_qbeg = _sa_qbeg;
            memcpy(sa_off, fields.a, sizeof(uint32_t) * 6);
            *sa_tid = bam_name2id(h, sa->s + alns.a[sa_idx]);
            *sa_pos = atoi(sa->s + alns.a[sa_idx] + sa_off[1]);
            *sa_is_rev = _sa_is_rev;
            *res = _res;
        }
    }
    free(alns.a);
    free(fields.a);
    free(sa_cigar.a);
    
    if (sa_idx >= 0) {
        return 1;
    } else {
        return 0;
    }
}
