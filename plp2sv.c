#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "plp2sv.h"
#include "cigar.h"
#include "array.h"

int plp2sv(bam_hdr_t *h, int tid, int pos, int n, int *n_plp, const bam_pileup1_t **plp, khash_t(sv_hash) *sv_h)
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
            int sa_is_rev, sa_tid, sa_pos, is_rev, qpos1, qpos2;
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

            qpos1 = is_rev ? cigar_res.clip[1] : cigar_res.clip[0];
            qpos1 += (sv1.ori1 ^ is_rev) ? cigar_res.qlen : 0;
            qpos2 = sa_is_rev ? sa_cig_res.clip[1] : sa_cig_res.clip[0];
            qpos2 += (sv1.ori2 ^ sa_is_rev) ? sa_cig_res.qlen : 0;
            //fprintf(stderr, "sa_qbeg = %d, sv_qbeg = %d, qpos1 = %d, qpos2 = %d\n", sa_qbeg, sv_qbeg, qpos1, qpos2);
            sv1.qdist = sa_qbeg < sv_qbeg ? qpos1 - qpos2 : qpos2 - qpos1;

            if (sv1.tid1 == sv1.tid2) {
                if (!(sa_is_rev ^ is_rev)) {
                    int _ori1 = sv1.ori1, _ori2 = sv1.ori2, rdist = sv1.pos2 - sv1.pos1;
                    if (sv1.pos1 > sv1.pos2) { _ori1 = sv1.ori2; _ori2 = sv1.ori1; rdist = sv1.pos1 - sv1.pos2; }
                    if (_ori1 & (!_ori2)) {
                        if (rdist > sv1.qdist) {
                            sv1.type = 'D';
                        } else {
                            sv1.type = 'I';
                        }
                    } else {
                        sv1.type = 'C';
                    }
                } else {
                    sv1.type = 'V';
                }
            } else {
                sv1.type = 'X';
            }

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
            } else {
                if (sv1.ori1 != kh_value(sv_h, k_iter).ori1 || sv1.ori2 != kh_value(sv_h, k_iter).ori2 ||
                        sv1.type != kh_value(sv_h, k_iter).type || sv1.qdist != kh_value(sv_h, k_iter).qdist) {
                    sv_t sv2 = kh_value(sv_h, k_iter);
                    fprintf(stderr, "SVs with matching breakpoints have differing values!\n");
                    fprintf(stderr, "SV1 ori1=%d, ori2=%d, type=%c, qdist=%d\n", sv1.ori1, sv1.ori2, sv1.type, sv1.qdist);
                    fprintf(stderr, "SV2 ori1=%d, ori2=%d, type=%c, qdist=%d\n", sv2.ori1, sv2.ori2, sv1.type, sv2.qdist);
                }
            }
        }
    }
    free(sup_alns.s);
    free(alns.a);
    free(sa_cigar.a);
    free(fields.a);
    return n_allele;
}
