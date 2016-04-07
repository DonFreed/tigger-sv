#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "plp2sv.h"
#include "cigar.h"
#include "array.h"

inline int plp2sv(bam_hdr_t *h, int tid, int pos, int n, int *n_plp, const bam_pileup1_t **plp, khash_t(sv_hash) *sv_h)
{
    int i, j, n_allele = 0, sv_qbeg, ret;
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
    khiter_t k_iter;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            int sa_is_rev, sa_tid, sa_pos, is_rev;
            b = plp[i][j].b;
            is_rev = !!(b->core.flag&BAM_FREVERSE);
            if (b->core.pos != pos && bam_endpos(b) != pos + 1) {
                continue; // Skip the middle of reads
            }
            is_front = b->core.pos == pos ? 1 : 0;
            tmp = bam_aux_get(b, "SA");
            if (!tmp) continue; // Skip reads without supplementary alignments
            n_cigar = b->core.n_cigar;
            parse_cigar(bam_get_cigar(b), n_cigar, &cigar_res);
            sv_qbeg = cigar_res.clip[is_rev];
            s = bam_aux2Z(tmp);
            sup_alns.l = 0;
            kputs(s, &sup_alns);
            ret = parse_sa_tag(h, &sup_alns, is_front, is_rev, sv_qbeg, &sa_cig_res, &sa_is_rev, &sa_qbeg, &sa_tid, &sa_pos);
            if (!ret) { continue; } // Supp alignments in wrong orientation
            sv1.tid1 = b->core.tid;
            sv1.tid2 = sa_tid;
            sv1.pos1 = pos;
            sv1.pos2 = (sa_qbeg < sv_qbeg) ^ sa_is_rev ? sa_pos + sa_cig_res.rlen - 2 : sa_pos -1;
            sv1.ori1 = !is_front;
            sv1.ori2 = (sa_qbeg < sv_qbeg) ^ sa_is_rev ? 1 : 0;
            /* FIXME: representing breakpoint positions as a 64-bit int may lead to collisions */
            if (sv1.pos1 < sv1.pos2) {
                sv1.id = (uint64_t)(sv1.tid1 & 0xff) << 56 | sv1.pos1 << 28 | sv1.pos2;
            } else {
                sv1.id = (uint64_t)(sv1.tid2 & 0xff) << 56 | sv1.pos2 << 28 | sv1.pos1;
            } 

            k_iter = kh_get(sv_hash, sv_h, sv1.id);
            if (k_iter == kh_end(sv_h)) {
                sv1.allele = n_allele++;
                k_iter = kh_put(sv_hash, sv_h, sv1.id, &ret);
                if (ret < 0) { return -1; }
                kh_value(sv_h, k_iter) = sv1;
            }
        }
    }
    free(sup_alns.s);
    free(alns.a);
    free(sa_cigar.a);
    free(fields.a);
    return n_allele;
}
