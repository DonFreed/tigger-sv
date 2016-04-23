#include <stdint.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "sv_qual.h"
#include "cigar.h"

int get_qual_data(bam_hdr_t *h, int tid, int pos, int n, int *n_plp,const bam_pileup1_t **plp, int n_alleles, khash_t(sv_hash) *sv_h, khash_t(sv_geno) *geno_h, mempool_t *mp)
{
    int i, j, nm, qbeg, res;
    bam1_t *b;
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
    qual_sum_t qual_sum;
    memset(&qual_sum, 0, sizeof(qual_sum));
    qual_sum.tid = tid;
    qual_sum.pos = pos;
    qual_sum.n_alleles = n_alleles;
    qual_sum.read_data = (uint16_t*)mp_alloc(mp, sizeof(uint16_t) * n * (n_alleles + 1));
    qual_sum.alleles = (allele_t*)mp_alloc(mp, sizeof(allele_t) * n_alleles);

    // read data //
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            int rd_idx, dp, allele, is_fwd;
            b = plp[i][j].b;
            qual_sum.mq_sum += b->core.qual * b->core.qual;
            parse_cigar(bam_get_cigar(b), b->core.n_cigar, &cigar_res);
            qual_sum.qlen_sum += cigar_res.qlen * cigar_res.qlen;
            dp = cigar_get_qual(bam_get_qual(b), pos - b->core.pos, bam_get_cigar(b), b->core.n_cigar);
            if (dp == UINT32_MAX) {
                fprintf(stderr, "Error parsing cigar, no quality score found\n");
                return -2;
            }
            nm = -1;
            if ((tmp = bam_aux_get(b, "NM"))) {
                float div;
                nm = bam_aux2i(tmp);
                div = (float)nm / (cigar_res.qlen + cigar_res.del);
                qual_sum.div_sum += div * div;
            } else {
                qual_sum.div_sum += (float)1.0;
            }
            if ((tmp = bam_aux_get(b, "AS"))) {
                int as = bam_aux2i(tmp);
                qual_sum.as_sum += as * as;
            } else {
                qual_sum.as_sum += 0;
            }
            if (b->core.pos != pos && bam_endpos(b) != pos + 1) {
                allele = 0;
                rd_idx = (i * n_alleles) + allele;
                qual_sum.read_data[rd_idx] += dp;
                continue;
            }
            tmp = bam_aux_get(b, "SA");
            if (!tmp) { 
                allele = 0; 
                rd_idx = (i * n_alleles) + allele;
                qual_sum.read_data[rd_idx] += dp;
                continue;
            }
            s = bam_aux2Z(tmp);
            sup_alns.l = 0;
            kputs(s, &sup_alns);
            is_fwd = !bam_is_rev(b);
            qbeg = cigar_res.clip[!is_fwd];
            res = parse_sa_tag(h, &sup_alns, b->core.pos == pos, !is_fwd, qbeg, &sa_cig_res, &sa_is_rev, &sa_qbeg, &sa_tid, &sa_pos);
            if (!res) { 
                allele = 0;
                rd_idx = (i * n_alleles) + allele;
                qual_sum.read_data[rd_idx] += dp;
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
                allele = kh_value(sv_h, k_iter).allele + 1;
                rd_idx = (i * n_alleles) + allele;
                qual_sum.read_data[rd_idx] += dp;
            } else {
                allele = 0;
                rd_idx = (i * n_alleles) + allele;
                qual_sum.read_data[rd_idx] += dp;
            }
        }
        qual_sum.n_reads += n_plp[i];
    }

    // allele data //
    for (k_iter = kh_begin(sv_h); k_iter != kh_end(sv_h); ++k_iter) {
        if (kh_exist(sv_h, k_iter)) {
            int allele = kh_value(sv_h, k_iter).allele;
            qual_sum.alleles[allele].tid = kh_value(sv_h, k_iter).tid2;
            qual_sum.alleles[allele].pos = kh_value(sv_h, k_iter).pos2;
            qual_sum.alleles[allele].ori1 = kh_value(sv_h, k_iter).ori1;
            qual_sum.alleles[allele].ori2 = kh_value(sv_h, k_iter).ori2;
            qual_sum.alleles[allele].type = kh_value(sv_h, k_iter).type;
            qual_sum.alleles[allele].qdist = kh_value(sv_h, k_iter).qdist;
        }
    }

    id = (uint64_t)tid << 32 | pos;
    k_iter = kh_put(sv_geno, geno_h, id, &res);
    if (res < 0) { return -1; }
    kh_value(geno_h, k_iter) = qual_sum;
            
    free(sup_alns.s);
    free(alns.a);
    free(fields.a);
    free(sa_cigar.a);
    return 1;
}

