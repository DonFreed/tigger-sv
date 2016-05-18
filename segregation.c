#include <stdio.h>
#include <math.h>
#include "segregation.h"

static double mendel_table[] = {
         //        child          fa  mo
         //  0/0   0/1   1/1
            1.00, 0.00, 0.00, // 0/0, 0/0
            0.75, 0.25, 0.00, // 0/0, 0/1
            0.00, 1.00, 0.00, // 0/0, 1/1
            0.75, 0.25, 0.00, // 0/1, 0/0
            0.25, 0.50, 0.25, // 0/1, 0/1
            0.00, 0.25, 0.75, // 0/1, 1/1
            0.00, 1.00, 0.00, // 1/1, 0/0
            0.00, 0.25, 0.75, // 1/1, 0/1
            0.00, 0.00, 1.00, // 1/1, 1/1
};

inline void update_likelihoods(double *ml, double *ol, int fa_gt, int mo_gt, int rr, int ra, int aa)
{
    // currently uses natural log
    double p_ref, p_het, p_alt;
    int m_rr, m_ra, m_aa;
    p_ref = mendel_table[fa_gt * 9 + mo_gt * 3];
    p_het = mendel_table[fa_gt * 9 + mo_gt * 3 + 1];
    p_alt = mendel_table[fa_gt * 9 + mo_gt * 3 + 2];
    *ol += lgamma((double)(rr + ra + aa + 1)) + 
            rr * log(p_ref) - lgamma(rr + 1) +
            ra * log(p_het) - lgamma(ra + 1) +
            aa * log(p_alt) - lgamma(aa + 1);
    m_rr = (int)(rr * p_ref + 0.499999);
    m_ra = (int)(ra * p_het + 0.499999);
    m_aa = (int)(aa * p_alt + 0.499999);
    *ml += lgamma((double)(rr + ra + aa + 1)) +
            m_rr * log(p_ref) - lgamma(m_rr + 1) +
            m_ra * log(p_het) - lgamma(m_ra + 1) +
            m_aa * log(p_alt) - lgamma(m_aa + 1);
}

double log_segregation(genotype_t *gts, int n, double mi_penalty, khash_t(ped) *ped_h)
{
    int i, j, trio_gts[27] = {0}; // 0-2 define child; 0,3,6 define mo, 0,9,18 define fa where 0=RR, 1=RA, 2=AA
    double ml, ol;
    khiter_t k;
    for (i = 0; i < n; ++i) {
        k = kh_get(ped, ped_h, i);
        if (k != kh_end(ped_h)) { // The sample is a child
            int idx = gts[i].gt;
            idx += (3 * gts[kh_value(ped_h, k).mo_col].gt);
            idx += (9 * gts[kh_value(ped_h, k).fa_col].gt);
            trio_gts[idx] += 1;
        }
    }
    ml = log(1.0); // maximum likelihood
    ol = log(1.0); // observed likelihood
    for (i = 0; i < 3; ++i) { // fa gt
        for (j = 0; j < 3; ++j) { // mo gt
            update_likelihoods(&ml, &ol, i, j, trio_gts[i * 9 + j * 3], trio_gts[i * 9 + j * 3 + 1], trio_gts[i * 9 + j * 3 + 2]);
            fprintf(stderr, "ol is now %f, ml is now %f\n", ol, ml);
        }
    }
    return ((ol - ml) / log(10)) * -10;
}
