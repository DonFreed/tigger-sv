#include <stdint.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "sv_qual.h"
#include "cigar.h"

inline int qual_push(qual_vec_t *quals, read_qual_t *qual1)
{
    if (quals->n >= quals->max) {
        quals->max = quals->n ? quals->n + 1 : 32;
        kroundup32(quals->max);
        if (!(quals->a = (read_qual_t*)realloc(quals->a, quals->max * sizeof(read_qual_t)))) {
            quals->a = 0; return 1;
        }
    }
    quals->a[quals->n++] = *qual1;
    return 0;
}

int get_qual_data(bam_hdr_t *h, int tid, int pos, int n, int *n_plp,const bam_pileup1_t **plp, khash_t(sv_hash) *sv_h, qual_vec_t *quals)
{
    int i, j, nm, qbeg, res;
    bam1_t *b;
    read_qual_t qual1;
    cigar_res_t cigar_res, sa_cig_res;
    uint8_t *tmp;
    const char *s;
    kstring_t sup_alns = {0, 0, 0};
    iarray_t alns = {0, 0, 0};
    iarray_t fields = {0, 0, 0};
    i32array_t sa_cigar = {0, 0, 0};
    int sa_is_rev, sa_qbeg, sa_tid, sa_pos, sa_bp_pos;
    uint64_t id;
    khiter_t k_iter;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            b = plp[i][j].b;
            qual1.qual = b->core.qual;
            parse_cigar(bam_get_qual(b), b->core.l_qseq, bam_get_cigar(b), b->core.n_cigar, &cigar_res);
            qual1.qlen = cigar_res.qlen;
            qual1.dp = cigar_get_qual(bam_get_qual(b), pos - b->core.pos, bam_get_cigar(b), b->core.n_cigar);
            if (qual1.dp == UINT32_MAX) {
                fprintf(stderr, "Error parsing cigar, no quality score found\n");
                return -2;
            }
            qual1.is_fwd = !bam_is_rev(b);
            nm = -1;
            if ((tmp = bam_aux_get(b, "NM"))) {
                nm = bam_aux2i(tmp);
                qual1.div = (float)nm / (qual1.qlen + cigar_res.del);
            } else {
                qual1.div = (float)1.0;
            }
            if (b->core.pos != pos && bam_endpos(b) != pos + 1) {
                qual1.allele = 0;
                if (qual_push(&quals[i], &qual1)) {
                    return -1;
                }
                continue;
            }
            tmp = bam_aux_get(b, "SA");
            if (!tmp) { 
                qual1.allele = 0; 
                if (qual_push(&quals[i], &qual1)) {
                    return -1;
                }
                continue;
            }
            s = bam_aux2Z(tmp);
            sup_alns.l = 0;
            kputs(s, &sup_alns);
            qbeg = cigar_res.clip[!qual1.is_fwd];
            res = parse_sa_tag(h, &sup_alns, b->core.pos == pos, !qual1.is_fwd, qbeg, &sa_cig_res, &sa_is_rev, &sa_qbeg, &sa_tid, &sa_pos);
            if (!res) { 
                qual1.allele = 0; 
                if (qual_push(&quals[i], &qual1)) {
                    return -1;
                }
                continue;
            }
            sa_bp_pos = (sa_qbeg < qbeg) ^ sa_is_rev ? sa_pos + sa_cig_res.rlen - 2 : sa_pos - 1;
            if (pos < sa_bp_pos) {
                id = (uint64_t)(b->core.tid & 0xff) << 56 | pos << 28 | sa_bp_pos;
            } else {
                id = (uint64_t)(sa_tid & 0xff) << 56 | sa_bp_pos << 28 | pos;
            }
            k_iter = kh_get(sv_hash, sv_h, id);
            if (k_iter != kh_end(sv_h)) {
                qual1.allele = kh_value(sv_h, k_iter).allele + 1;
                if (qual_push(&quals[i], &qual1)) {
                    return -1;
                }
            } else {
                qual1.allele = 0;
                if (qual_push(&quals[i], &qual1)) {
                    return -1;
                }
            }
        }
    }
    free(sup_alns.s);
    free(alns.a);
    free(fields.a);
    free(sa_cigar.a);
    return 1;
}

