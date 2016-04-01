#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "plp2sv.h"

inline void parse_cigar(uint8_t *qual, int32_t l_qseq, const uint32_t *cigar, int n_cigar, int clip[2], int *qlen, int *rlen, int *ins, int *del, int clip_q[2])
{
    int k;
    if (qual) {
        clip_q[0] = qual[0]; 
        clip_q[1] = qual[l_qseq];
    }
    clip[0] = 0; clip[1] = 0; *qlen = 0; *rlen = 0; *ins = 0; *del = 0;
    for (k = 0; k < n_cigar; ++k) {
        int op = bam_cigar_op(cigar[k]);
        int oplen = bam_cigar_oplen(cigar[k]);
        if (oplen == 0) continue;
        if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP) *qlen += oplen;
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            clip[!!k] = oplen;
            if (qual && op == BAM_CSOFT_CLIP) {
                clip_q[!!k] = !k? qual[oplen] : qual[clip[0] + *qlen - 1];
            }
        }
        if (bam_cigar_type(op)&2) *rlen += oplen;
        if (op == BAM_CINS) *ins += oplen;
        else if (op == BAM_CDEL) *del += oplen;
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

inline int plp2sv(int tid, int pos, int n, int *n_plp, const bam_pileup1_t **plp, sv_vec_t *sv)
{
    int i, j, k;
    uint8_t *tmp;
    bam1_t *b;
    const char *s;
    int sa_clip[2], sa_qbeg, sa_qlen, sa_rlen, is_front, n_cigar;
    kstring_t sup_alns = {0, 0, 0}; // FIXME: frequent memory allocation
    iarray_t alns = {0, 0, 0}; // FIXME: frequent memory allocation
    iarray_t fields = {0, 0, 0};
    i32array_t sa_cigar = {0, 0, 0}; // FIXME: frequent memory allocation
    sv_t sv1;
    sv->n = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            int sa_idx = -1, lqbeg = -1, sa_off[8], sa_ins, sa_del, clip_q[2], sa_clip_q[2], sa_n_cigar;
            int is_rev;
            b = plp[i][j].b;
            is_rev = !!(b->core.flag&BAM_FREVERSE);
            if (b->core.pos != pos && bam_endpos(b) != pos + 1) {
                continue; // Skip the middle of reads
            }
            is_front = b->core.pos == pos ? 1 : 0;
            tmp = bam_aux_get(b, "SA");
            if (!tmp) continue; // Skip reads without supplementary alignments
            n_cigar = b->core.n_cigar;
            parse_cigar(bam_get_qual(b), b->core.l_qseq, bam_get_cigar(b), n_cigar, sv1.clip, &sv1.qlen, &sv1.rlen, &sv1.ins, &sv1.del, clip_q);
            sv1.qbeg = sv1.clip[is_rev];
            s = bam_aux2Z(tmp);
            sup_alns.l = 0;
            kputs(s, &sup_alns);
            alns.n = ksplit_core(sup_alns.s, ';', &alns.max, &alns.a);
            for (k = 0; k < alns.n; ++k) {
                fields.n = ksplit_core(sup_alns.s + alns.a[k], ',', &fields.max, &fields.a);
                sa_cigar.n = 0;
                str2cigar(sup_alns.s + alns.a[k] + fields.a[3], &sa_cigar, &sa_n_cigar);
                parse_cigar(0, 0, sa_cigar.a, sa_n_cigar, sa_clip, &sa_qlen, &sa_rlen, &sa_ins, &sa_del, sa_clip_q);
                fprintf(stderr, "Supp alignment has orientation %c with sa_clip %d, %d\n", *(sup_alns.s + alns.a[k] + fields.a[2]), sa_clip[0], sa_clip[1]);
                sa_qbeg = *(sup_alns.s + alns.a[k] + fields.a[2]) == '-' ? sa_clip[1] : sa_clip[0];
                // Is the supp alignment in the correct orientation?
                fprintf(stderr, "qbeg = %d, sa_qbeg = %d, is_front = %d, is_rev = %d\n", sv1.qbeg, sa_qbeg, is_front, is_rev);
                if (is_front ^ (!is_rev)) {
                    if (sv1.qbeg < sa_qbeg) {
                        if (sa_idx >= 0) {
                            if (sa_qbeg < lqbeg) {
                                sa_idx = k;
                                lqbeg = sa_qbeg;
                                memcpy(&sa_off, &fields.a, sizeof(int) * 6);
                            }
                        } else {
                            sa_idx = k;
                            lqbeg = sa_qbeg;
                            memcpy(&sa_off, &fields.a, sizeof(int) * 6);
                        }
                    }
                } else {
                    if (sv1.qbeg > sa_qbeg) {
                        if (sa_idx >= 0) {
                            if (sa_qbeg > lqbeg) {
                                sa_idx = k;
                                lqbeg = sa_qbeg;
                                memcpy(&sa_off, &fields.a, sizeof(int) * 6);
                            }
                        } else {
                            sa_idx = k;
                            lqbeg = sa_qbeg;
                            memcpy(&sa_off, &fields.a, sizeof(int) * 6);
                        }
                    }
                }
            }
            if (sa_idx < 0) { fprintf(stderr, "Wrong orientation\n"); continue; } // Supp alignments in wrong orientation
            sv1.to_pos = (uint32_t)atoi(sup_alns.s + alns.a[sa_idx] + sa_off[1]);
            sv1.tid = b->core.tid;
            sv1.pos = b->core.pos;
            sv1.id = sv1.pos < sv1.to_pos ? (uint64_t)sv1.pos << 32 | sv1.to_pos : (uint64_t)sv1.to_pos << 32 | sv1.pos;
            sv1.flag = b->core.flag;
            sv1.mapq = b->core.qual;
            sv1.nm = -1; sv1.diff = -1; sv1.score = -1;
            if ((tmp = bam_aux_get(b, "NM"))) {
                sv1.nm = bam_aux2i(tmp);
                sv1.diff = (double)sv1.nm / (sv1.qlen + sv1.del);
            }
            if ((tmp =  bam_aux_get(b, "AS"))) {
                sv1.score = bam_aux2i(tmp);
            }
            sv1.tipq = is_front? clip_q[0] : clip_q[1];
            if (sv->n >= sv->max) {
                sv->max = sv->n ? sv->n + 1 : 16;
                kroundup32(sv->max);
                sv->sv = (sv_t*)realloc(sv->sv, sv->max * sizeof(sv_t));
               if (!sv->sv) return -1;
            }
            sv->sv[sv->n++] = sv1;
        }
    }
    free(sup_alns.s);
    free(alns.a);
    free(sa_cigar.a);
    free(fields.a);
    return sv->n;
}
