#include <stdio.h>
#include <math.h>
#include "segregation.h"

void update_likelihoods(double *ml, double *ol, int fa_gt, int mo_gt, int rr, int ra, int aa);
{
    // currently uses natural log
    double _ml, tmp;
    int i, j, m_rr, m_ra, m_aa;
    *ol += lgamma((double)(rr + ra + aa + 1)) + 
            rr * log(p) - lgamma(rr + 1) +
            ra * log(h) - lgamma(ra + 1) +
            aa * log(q) - lgamma(aa + 1);
    m_rr = (int)(rr * p + 0.499999);
    m_ra = (int)(ra * h + 0.499999);
    m_aa = (int)(aa * q + 0.499999);
    *ml += lgamma((double)(rr + ra + aa + 1)) +
            m_rr * log(p) - lgamma(m_rr + 1) +
            m_ra * log(h) - lgamma(m_ra + 1) +
            m_aa * log(q) - lgamma(m_aa + 1);
}

double log_segregation(genotype_t *gts, int n, double mi_penalty, ped_t *ped)
{
    int i, j, trio_gts[27]; // 0-2 define child; 0,3,6 define mo, 0,9,18 define fa where 0=RR, 1=RA, 2=AA
    double ml, ol;
    for (i = 0; i < n; ++i) {
        if (sample[i] is child) {
            int idx = gts[i].gt;
            idx += (3 * gts[mo_idx]);
            idx += (9 * gts[fa_idx]);
            trio_gts[idx] += 1;
        }
    }
    ml = log10(1.0); // maximum likelihood
    ol = lot10(1.0); // observed likelihood
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            update_likelihoods(&ml, &ol, i, j, trio_gts[i * 9 + i * 3], trio_gts[i * 9 + i * 3 + 1], trio_gts[i * 9 + i * 3 + 2]);
        }
    }
    return ((ol - ml) / log(10)) * -10;
}
