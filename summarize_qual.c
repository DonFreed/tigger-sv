#include <string.h>
#include "summarize_qual.h"

int summarize_qual(uint32_t tid, uint32_t pos, int n, int n_alleles, khash_t(sv_hash) *sv_h, qual_vec_t *quals, mempool_t *mp, khash_t(sv_geno) *geno_h)
{
    qual_sum_t qual_sum;
    int i, j, ret;
    khiter_t k_iter;
    sv_t sv1;
    uint64_t id = (uint64_t)tid << 32 | pos;
    memset(&qual_sum, 0, sizeof(qual_sum));
    qual_sum.tid = tid;
    qual_sum.pos = pos;
    qual_sum.n_alleles = n_alleles;
    qual_sum.read_data = (uint16_t*)mp_alloc(mp, sizeof(uint16_t) * n * (n_alleles + 1));
    qual_sum.alleles = (allele_t*)mp_alloc(mp, sizeof(allele_t) * n_alleles);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < quals[i].n; ++j) {
            int rd_idx = (i * n_alleles) + quals[i].a[j].allele;
            qual_sum.qlen_sum += quals[i].a[j].qlen;
            qual_sum.div_sum += quals[i].a[j].div;
            qual_sum.mq_sum += quals[i].a[j].qual;
            qual_sum.read_data[rd_idx] += quals[i].a[j].dp;
        }
        qual_sum.n_reads += quals[i].n;
    }
    for (k_iter = kh_begin(sv_h); k_iter != kh_end(sv_h); ++k_iter) {
        if (kh_exist(sv_h, k_iter)) {
            int allele;
            sv1 = kh_value(sv_h, k_iter);
            allele = sv1.allele;
            qual_sum.alleles[allele].tid = sv1.tid2;
            qual_sum.alleles[allele].pos = sv1.pos2;
            qual_sum.alleles[allele].ori1 = sv1.ori1;
            qual_sum.alleles[allele].ori2 = sv1.ori2;
        }
    }
    k_iter = kh_put(sv_geno, geno_h, id, &ret);
    if (ret < 0) { return -1; }
    kh_value(geno_h, k_iter) = qual_sum;
    return 1;
}

