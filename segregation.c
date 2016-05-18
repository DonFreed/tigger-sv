#include <stdio.h>
#include <math.h>
#include "segregation.h"

void gt_to_alleles(int gt, int *alleles)
{
    int i;
    for (i = 0;; ++i) {
        if ( (i + 1) * (i + 2) / 2 > gt) {
            break;
        }
    }
    alleles[0] = i - (( (i + 1) * (i + 2) / 2 ) - gt - 1);
    alleles[1] = i;
    return;
}

inline void mendel_probs(int fa, int mo, double *probs, double mi_prob)
{
    int i, j;
    char fa_a[2], mo_a[2];
    memset(probs, 0, sizeof(double) * 3);
    fa_a[0] = fa & 2; fa_a[1] = (fa & 1 || fa & 2);
    mo_a[0] = mo & 2; mo_a[1] = (mo & 1 || mo & 2);
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            probs[fa_a[i] + mo_a[j]] += 0.25;
        }
    }
    for (i = 0; i < 3; ++i) {
        if (probs[i] < 0.1) { // will cause problems with mi_prob > 0.1
            if (i == 0) {
                if (probs[1] > 0.1) {
                    probs[0] += mi_prob;
                    probs[1] -= mi_prob / 2;
                    probs[2] -= mi_prob / 2;
                } else {
                    probs[0] += mi_prob / 2;
                    probs[2] -= mi_prob / 2;
                }
            } else if (i == 1) {
                if (probs[0] > 0.1) {
                    probs[1] += mi_prob;
                    probs[0] -= mi_prob;
                } else {
                    probs[1] += mi_prob;
                    probs[2] -= mi_prob;
                }
            } else {
                if (probs[1] > 0.1) {
                    probs[2] += mi_prob;
                    probs[0] -= mi_prob / 2;
                    probs[1] -= mi_prob / 2;
                } else {
                    probs[2] += mi_prob / 2;
                    probs[0] -= mi_prob / 2;
                }
            }
        }
    }
}

inline void update_likelihoods(double *ml, double *ol, double *probs, int rr, int ra, int aa)
{
    int m_rr, m_ra, m_aa;
    *ol += lgamma((double)(rr + ra + aa + 1)) + 
            rr * log(probs[0]) - lgamma(rr + 1) +
            ra * log(probs[1]) - lgamma(ra + 1) +
            aa * log(probs[2]) - lgamma(aa + 1);
    m_rr = (int)(rr * probs[0] + 0.499999);
    m_ra = (int)(ra * probs[1] + 0.499999);
    m_aa = (int)(aa * probs[2] + 0.499999);
    *ml += lgamma((double)(rr + ra + aa + 1)) +
            m_rr * log(probs[0]) - lgamma(m_rr + 1) +
            m_ra * log(probs[1]) - lgamma(m_ra + 1) +
            m_aa * log(probs[2]) - lgamma(m_aa + 1);
}

double log_segregation(genotype_t *gts, int n, double mi_prob, khash_t(ped) *ped_h)
{
    int i, j, trio_gts[27] = {0}; // 0-2 define child; 0,3,6 define mo, 0,9,18 define fa where 0=RR, 1=RA, 2=AA
    double ml, ol, probs[3];
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
            mendel_probs(i, j, probs, mi_prob);
            update_likelihoods(&ml, &ol, probs, trio_gts[i * 9 + j * 3], trio_gts[i * 9 + j * 3 + 1], trio_gts[i * 9 + j * 3 + 2]);
            fprintf(stderr, "ol is now %f, ml is now %f\n", ol, ml);
        }
    }
    return ((ol - ml) / log(10)) * -10;
}

int n_mi(genotype_t *gts, int n, khash_t(ped) *ped_h)
{
    int i, j, fa_a[2], mo_a[2], child_a[2], ret = 0;
    khiter_t k;
    for (i = 0; i < n; ++i) {
        k = kh_get(ped, ped_h, i);
        if (k != kh_end(ped_h)) { // Father and mother are present
            gt_to_alleles(gts[kh_value(ped_h, k).fa_col].gt, fa_a);
            gt_to_alleles(gts[kh_value(ped_h, k).mo_col].gt, mo_a);
            gt_to_alleles(gts[i].gt, child_a);
            for (j = 0; j < 2; ++j) {
                int mc = 0; // mendelian consistency
                for (k = 0; k < 2; ++k) {
                    if (child_a[j] == fa_a[k] || child_a[j] == mo_a[k]) {
                        mc = 1;
                    }
                }
                if (!mc) {
                    ret += 1;
                }
            }
        }
    }
    return ret;
}

