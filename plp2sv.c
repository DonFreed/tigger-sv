#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"

inline void parse_cigar(const uint32_t *cigar, int *clip[2], int *qbeg, int *qlen, int *rlen, int *ins, int *del, int *clip_q[2])
{
    pass;
}

inline void str2cigar(const char *s, i32array_t *cigar)
{
    pass;
}

inline int plp2sv(int tid, int pos, int n; int *n_plp; bam_pileup1_t *plp; sv_vec_t *sv)
{
    int i, j, k;
    uint8_t *tmp;
    bam1_t *b;
    const char *s;
    int clip[2], qbeg, qlen, rlen, sa_clip[2], sa_qbeg, sa_qlen, sa_rlen, is_front;
    int fields_n, fields_m = 8, fields_off[8];
    kstring_t sup_alns = {0, 0, 0}; // FIXME: frequent memory allocation
    iarray_t alns = {0, 0, 0}; // FIXME: frequent memory allocation
    i32array_t sa_cigar = {0, 0, 0}; // FIXME: frequent memory allocation
    sv_t sv1;
    sv->n = 0;
    for (i = 0; i < ni; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            int sa_idx = -1, lqbeg = -1, sa_off[8], sa_ins, sa_del, clip_q[2], sa_clip_q[2];
            int is_rev;
            b = plp[i][j].b;
            is_rev = !!(b->core.flag&BAM_FREVERSE)
            if (b->core.pos != pos && bam_endpos(b) != pos + 1) {
                continue; // Skip the middle of reads
            }
            is_front = b->core.pos == pos ? 1 : 0;
            tmp = bam_aux_get(b, "SA");
            if (!tmp) continue; // Skip reads without supplementary alignments
            parse_cigar(bam_get_cigar(b), &sv1.clip, &sv1.qbeg, &sv1.qlen, &sv1.rlen, &sv1.ins, &sv1.del, &clip_q);
            s = bam_aux2Z(tmp);
            sup_alns.n = 0;
            kputs(s, &sup_alns);
            alns.n = ksplit_core(sup_alns.s, ';', &aln.max, &alns.a);
            for (k = 0; k < alns.n; ++k) {
                fields.n = ksplit_core(sup_alns.s + alns.offsets[k], ',', fields_m, &fields_off);
                str2cigar(sup_alns.s + alns.offsets[k] + fields_off[3], sa_cigar);
                parse_cigar(sa_cigar.a, &sa_clip, &sa_qbeg, &sa_qlen, &sa_rlen, &sa_ins, &sa_del, &sa_clip_q);
                // Is the supp alignment in the correct orientation?
                if (is_front ^ is_rev) {
                    if (sv1.qbeg < sa_qbeg) {
                        if (sa_idx >= 0) {
                            if (sa_qbeg < lqbeg) {
                                sa_idx = k;
                                lqbeg = sa_qbeg;
                                memcpy(&sa_off, &fields_off, sizeof(fields_off));
                            }
                        } else {
                            sa_idx = k;
                            lqbeg = sa_qbeg;
                            memcpy(&sa_off, &fields_off, sizeof(fields_off));
                        }
                    }
                } else {
                    if (sv1.qbeg > sa_qbeg) {
                        if (sa_idx >= 0) {
                            if (sa_qbeg > lqbeg) {
                                sa_idx = k;
                                lqbeg = sa_qbeg;
                                memcpy(&sa_off, &fields_off, sizeof(fields_off));
                            }
                        } else {
                            sa_idx = k;
                            lqbeg = sa_qbeg;
                            memcpy(&sa_off, &fields_off, sizeof(fields_off));
                        }
                    }
                }
            }
            if (sa_idx < 0) continue; // Supp alignments in wrong orientation
            qual = bam_get_qual(b);
            sv1.to_pos = (uint32_t)atoi(sup_alns.s + alns.offsets[sa_idx] + sa_off[1]);
            sv1.tid = b->core.tid;
            sv1.pos = b->core.pos;
            sv1.id = sv1.pos < sv1.to_pos ? sv1.pos << 32 | sv1.to_pos : sv1.to_pos << 32 | sv1.pos;
            sv1.flag = b->core.flag;
            sv1.mapq = b->core.qual;
            sv1.nm = -1; sv1.diff = -1; sv1.score = -1;
            if (tmp = bam_aux_get(b, "NM")) {
                sv1.nm = bam_aux2i(tmp);
                sv1.diff = (double)sv1.nm / (sv1.qlen + sv1.del);
            }
            if (tmp =  bam_aux_get(b, "AS")) {
                sv1.score = bam_aux2i(tmp);
            }
            sv1.tipq = is_front? clip_q[0] : clip_q[1]
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
    return sv->n;
}
