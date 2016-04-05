#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "plp2sv.h"
#include "cigar.h"
#include "array.h"

inline int plp2sv(bam_hdr_t *h, int tid, int pos, int n, int *n_plp, const bam_pileup1_t **plp, sv_vec_t *sv)
{
    int i, j, k, l, sv_qbeg;
    uint8_t *tmp;
    bam1_t *b;
    const char *s;
    int sa_qbeg, is_front, n_cigar;
    kstring_t sup_alns = {0, 0, 0}; // FIXME: frequent memory allocation
    iarray_t alns = {0, 0, 0}; // FIXME: frequent memory allocation
    iarray_t fields = {0, 0, 0};
    i32array_t sa_cigar = {0, 0, 0}; // FIXME: frequent memory allocation
    sv_t sv1;
    cigar_res_t cigar_res, sa_cig_res;
    sv->n = 0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            uint32_t sa_off[8] = {0};
            int sa_idx = -1, sa_is_rev, lqbeg = -1, sa_n_cigar, lpos, lrlen, lis_rev;
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
            parse_cigar(bam_get_qual(b), b->core.l_qseq, bam_get_cigar(b), n_cigar, &cigar_res);
            sv_qbeg = cigar_res.clip[is_rev];
            s = bam_aux2Z(tmp);
            sup_alns.l = 0;
            kputs(s, &sup_alns);
            alns.n = ksplit_core(sup_alns.s, ';', &alns.max, &alns.a);
            for (k = 0; k < alns.n; ++k) {
                fields.n = ksplit_core(sup_alns.s + alns.a[k], ',', &fields.max, &fields.a);
                sa_cigar.n = 0;
                str2cigar(sup_alns.s + alns.a[k] + fields.a[3], &sa_cigar, &sa_n_cigar);
                parse_cigar(0, 0, sa_cigar.a, sa_n_cigar, &sa_cig_res);
                fprintf(stderr, "Supp alignment has orientation %c with sa_clip %d, %d\n", *(sup_alns.s + alns.a[k] + fields.a[2]), sa_cig_res.clip[0], sa_cig_res.clip[1]);
                sa_is_rev = *(sup_alns.s + alns.a[k] + fields.a[2]) == '-';
                sa_qbeg = sa_is_rev ? sa_cig_res.clip[1] : sa_cig_res.clip[0];
                // Is the supp alignment in the correct orientation?
                fprintf(stderr, "qbeg = %d, sa_qbeg = %d, is_front = %d, is_rev = %d\n", sv_qbeg, sa_qbeg, is_front, is_rev);
                if (is_front ^ (!is_rev)) {
                    if (sv_qbeg < sa_qbeg) {
                        if (sa_idx >= 0) {
                            if (sa_qbeg < lqbeg) {
                                sa_idx = k;
                                lqbeg = sa_qbeg;
                                memcpy(sa_off, fields.a, sizeof(uint32_t) * 6);
                                lpos = atoi(sup_alns.s + alns.a[sa_idx] + sa_off[1]); 
                                lrlen = sa_cig_res.rlen -1;
                                lis_rev = sa_is_rev;
                            }
                        } else {
                            sa_idx = k;
                            lqbeg = sa_qbeg;
                            memcpy(sa_off, fields.a, sizeof(uint32_t) * 6);
                            lpos = atoi(sup_alns.s + alns.a[sa_idx] + sa_off[1]); 
                            lrlen = sa_cig_res.rlen -1;
                            lis_rev = sa_is_rev;
                        }
                    }
                } else {
                    if (sv_qbeg > sa_qbeg) {
                        if (sa_idx >= 0) {
                            if (sa_qbeg > lqbeg) {
                                sa_idx = k;
                                lqbeg = sa_qbeg;
                                memcpy(sa_off, fields.a, sizeof(uint32_t) * 6);
                                lpos = atoi(sup_alns.s + alns.a[sa_idx] + sa_off[1]);
                                lrlen = sa_cig_res.rlen -1;
                                lis_rev = sa_is_rev;
                            }
                        } else {
                            sa_idx = k;
                            lqbeg = sa_qbeg;
                            memcpy(sa_off, fields.a, sizeof(uint32_t) * 6);
                            lpos = atoi(sup_alns.s + alns.a[sa_idx] + sa_off[1]);
                            lrlen = sa_cig_res.rlen -1;
                            lis_rev = sa_is_rev;
                        }
                    }
                }
            }
            if (sa_idx < 0) { fprintf(stderr, "Wrong orientation\n"); continue; } // Supp alignments in wrong orientation
            sv1.tid1 = b->core.tid;
            sv1.tid2 = bam_name2id(h, sup_alns.s + alns.a[sa_idx]);
            sv1.pos1 = pos;
            sv1.pos2 = (lqbeg < sv_qbeg) ^ lis_rev ? lpos + lrlen - 1 : lpos;
            sv1.ori1 = !is_front;
            sv1.ori2 = (lqbeg < sv_qbeg) ^ lis_rev ? 1 : 0;
            sv1.id = sv1.pos1 < sv1.pos2 ? (uint64_t)sv1.pos1 << 32 | sv1.pos2 : (uint64_t)sv1.pos2 << 32 | sv1.pos1;
            
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