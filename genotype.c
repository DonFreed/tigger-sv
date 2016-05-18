#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include "tigger_version.h"
#include "genotype.h"
#include "asa147.h"
#include "segregation.h"

genotype_t *gts_init(int n, int n_var)
{
    genotype_t *gts = (genotype_t*)calloc(n, sizeof(genotype_t));
    int i, n_pl = (n_var + 1) * (n_var + 2) / 2;
    int *pl = (int*)calloc(n * n_pl, sizeof(int));
    fprintf(stderr, "Initalizing\n");
    for (i = 0; i < n; ++i, pl += n_pl) {
        gts[i].pl = pl;
    }
    return gts;
}

void gts_destroy(genotype_t *gts)
{
    free(gts[0].pl);
    free(gts);
}

void print_header(bam_hdr_t *h, int optind, int n, char *argv[])
{
    int i;
    printf("##fileformat=VCFv4.2\n");
    printf("##source=tigger-sv-%s\n", TIGGERSV_VERSION);
    //printf("##ALT=<ID=DEL,Description=\"Deletion relative to the reference genome\">\n");
    //printf("##ALT=<ID=INS,Description=\"Insertion relative to the reference genome\">\n");
    //printf("##ALT=<ID=INV,Description=\"Inversion relative to the reference genome\">\n");
    //printf("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the structural variant along the reference\">\n");
    //printf("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
    printf("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of the variant\">\n");
    printf("##INFO=<ID=QGAP,Number=A,Type=Integer,Description=\"Length of the gap on the query sequence\">\n");
    printf("##INFO=<ID=MAPQ,Number=1,Type=Float,Description=\"RMS mapping quality of supporing alignments\">\n");
    printf("##INFO=<ID=DIV,Number=1,Type=Float,Description=\"RMS divergence of supporting alignments\">\n");
    printf("##INFO=<ID=QLEN,Number=1,Type=Float,Description=\"RMS length of supporting alignments along the query\">\n");
    printf("##INFO=<ID=SC,Number=1,Type=Integer,Description=\"RMS alignment score of supporting alignments\">\n");
    printf("##INFO=<ID=AS_SC,Number=R,Type=Float,Description=\"RMS alignment score of alignments supporting the alternate allele(s)\"\n");
    //printf("##INFO=<ID=TIPQ,Number=1,Type=Integer,Description=\"RMS quality/depth of the supporting alignments\">\n");
    printf("##INFO=<ID=BPD,Number=.,Type=Integer,Description=\"Depths at breakpoint(s). All reference alleles are listed first followed by alternate allele pairs\">\n");
    printf("##FORMAT=<ID=BPD,Number=.,Type=Integer,Description=\"Depths at breakpoint(s). All reference alleles are listed first followed by alternate allele pairs\">\n");
    printf("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n");
    printf("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the nearest integer\">\n");
    for (i = 0; i < h->n_targets; ++i) {
        printf("##contig=<ID=%s,length=%" PRIu32 ">\n", h->target_name[i], h->target_len[i]);
    }
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (i = 0; i < n; ++i) {
        printf("\t%s", argv[optind + i]);
    }
    putchar('\n');
}

double p_err = 0.05, p_het = 0.5, p_hom = 0.99;

inline void genotype_samples(genotype_t *gts, qual_sum_t *qual1, qual_sum_t **qual2, int *var_idx, int n_var, int n) {
    int n_pl = (n_var + 1) * (n_var + 2) / 2;
    int i, j, k, l, min_pl, pl_idx, next_min_pl;
    uint16_t *dp1, *dp2;
    double *pl = (double*)malloc(sizeof(double) * n_pl);
    for (l = 0; l < n; ++l) {
        memset(pl, 0.0, sizeof(double) * n_pl);
        dp1 = qual1->read_data + qual1->n_alleles * l;
        for (pl_idx = 0, i = 0; i <= n_var; ++i) {
            for (j = 0; j <= i; ++j, ++pl_idx) {
                if (i == 0 && j == 0) {
                    pl[pl_idx] += dp1[0] * log10(1.0 - p_err);
                    for (k = 0; k < n_var; ++k) {
                        dp2 = qual2[k]->read_data + qual2[k]->n_alleles * l;
                        pl[pl_idx] += dp2[0] * log10(1.0 - p_err);
                        pl[pl_idx] += dp1[var_idx[k * 2] + 1] * log10(p_err);
                        pl[pl_idx] += dp2[var_idx[k * 2 + 1] + 1] * log10(p_err);
                    }
                } else if (i == j) {
                    pl[pl_idx] += dp1[0] * log10(1.0 - p_hom);
                    for (k = 0; k < n_var; ++k) {
                        dp2 = qual2[k]->read_data + qual2[k]->n_alleles * l;
                        pl[pl_idx] += dp2[0] * log10(1.0 - p_hom);
                        if (k + 1 == i) {
                            pl[pl_idx] += dp1[var_idx[k * 2] + 1] * log10(p_hom);
                            pl[pl_idx] += dp2[var_idx[k * 2 + 1] + 1] * log10(p_hom);
                        } else {
                            pl[pl_idx] += dp1[var_idx[k * 2] + 1] * log10(p_err);
                            pl[pl_idx] += dp2[var_idx[k * 2 + 1] + 1] * log10(p_err);
                        }
                    }
                } else {
                    if (j == 0) {
                        pl[pl_idx] += dp1[0] * log10(1.0 - p_het);
                    } else {
                        pl[pl_idx] += dp1[0] * log10(1.0 - p_hom);
                    }
                    for (k = 0; k < n_var; ++k) {
                        dp2 = qual2[k]->read_data + qual2[k]->n_alleles * l;
                        if (j == 0) {
                            pl[pl_idx] += dp2[0] * log10(1.0 - p_het);
                        } else {
                            pl[pl_idx] += dp2[0] * log10(1.0 - p_hom);
                        }
                        if (k + 1 == i || k + 1 == j) {
                            pl[pl_idx] += dp1[var_idx[k * 2] + 1] * log10(p_het);
                            pl[pl_idx] += dp2[var_idx[k * 2 + 1] + 1] * log10(p_het);
                        } else {
                            pl[pl_idx] += dp1[var_idx[k * 2] + 1] * log10(p_err);
                            pl[pl_idx] += dp2[var_idx[k * 2 + 1] + 1] * log10(p_err);
                        }
                    }
                }
                fprintf(stderr, "i is %d, j is %d, pl_idx is %d, pl=[%f,%f,%f]\n", i, j, pl_idx, pl[0], pl[1], pl[2]);
            }
        }
        for (i = 0; i < n_pl; ++i) {
            pl[i] = pl[i] * -10;
        }
        for (i = 1, pl_idx = 0, min_pl = pl[0], next_min_pl = INT_MAX; i < n_pl; ++i) {
            if (pl[i] < min_pl) {
                pl_idx = i;
                next_min_pl = min_pl;
                min_pl = pl[i];
            } else if (pl[i] < next_min_pl) {
                next_min_pl = pl[i];
            }
        }
        gts[l].gt = pl_idx;
        gts[l].genotype_confidence = (int)next_min_pl - (int)min_pl;
        for (i = 0; i < n_pl; ++i) {
            gts[l].pl[i] = (int)pl[i] - (int)min_pl;
        }
    }
    free(pl);
    return;
}

inline void print_genotypes(genotype_t *gts, qual_sum_t *qual1, qual_sum_t **qual2, int *var_idx, int n_var, int n) {
    int n_pl = (n_var + 1) * (n_var + 2) / 2;
    int i, j, k, l;
    uint16_t *dp1, *dp2;
    printf("\tGT:BPD:GQ:PL");
    for (l = 0; l < n; ++l) {
        putchar('\t');
        dp1 = qual1->read_data + qual1->n_alleles * l;
        for (i = 0;; ++i) {
            if ( (i + 1) * (i + 2) / 2 > gts[l].gt) {
                fprintf(stderr, "with i=%d, %d is greater than %d\n", i, (i + 1) * (i + 2) / 2, gts[l].gt);
                break;
            }
        }
        j = i - (( (i + 1) * (i + 2) / 2 ) - gts[l].gt - 1);
        printf("%d/%d:", j, i);

        printf("%d", dp1[0]);
        for (k = 0; k < n_var; ++k) {
            dp2 = qual2[k]->read_data + qual2[k]->n_alleles * l;
            printf(",%d", dp2[0]);
        }
        for (k = 0; k < n_var; ++k) {
            dp2 = qual2[k]->read_data + qual2[k]->n_alleles * l;
            printf(",%d,%d", dp1[var_idx[k * 2] + 1], dp2[var_idx[k * 2 + 1] + 1]);
        }
        printf(":%d:", gts[l].genotype_confidence);
        for (i = 0; i < n_pl; ++i) {
            if (i) putchar(',');
            printf("%d", gts[l].pl[i]);
        }
    }
    return;
}

int genotype_sv(bam_hdr_t *h, int n, khash_t(sv_geno) *geno_h, int min_dp, khash_t(ped) *ped_h, double mi_prob)
{
    int i, j, l;
    khiter_t k1, k2;
    qual_sum_t *qual1, *qtmp, **qual2 = 0;
    allele_t *a1, *a2;
    int n_var, m_var = 0;
    int *var_idx = 0, allele_dp;
    uint16_t *dp1, *dp2;
    uint64_t id;
    genotype_t *gts;

    for (k1 = kh_begin(geno_h); k1 != kh_end(geno_h); ++k1) {
        if (kh_exist(geno_h, k1)) {
            qual1 = &(kh_value(geno_h, k1));
            a1 = qual1->alleles;
            for (i = 0, n_var = 0; i + 1 < qual1->n_alleles; ++i) {
                if (a1[i].genotyped) { fprintf(stderr, "Allele %d genotyped\n", i); continue; }
                id = (uint64_t)a1[i].tid << 32 | a1[i].pos;
                k2 = kh_get(sv_geno, geno_h, id);
                if (k2 == kh_end(geno_h)) {
                    fprintf(stderr, "Error: breakpoint mate of %d:%d not found at %d:%d.\n", (int)qual1->tid, (int)qual1->pos, (int)a1[i].tid, (int)a1[i].pos);
                    continue;
                }
                qtmp = &(kh_value(geno_h, k2));
                a2 = qtmp->alleles;
                for (j = 0; j + 1 < qtmp->n_alleles; ++j) {
                    if (qual1->tid == a2[j].tid && qual1->pos == a2[j].pos) break;
                }
                if (j + 1 == qtmp->n_alleles) {
                    fprintf(stderr, "Error: breakpoint does not have a corresponding allele from %d:%d to %d:%d.\n", (int)qual1->tid, (int)qual1->pos, (int)a1[i].tid, (int)a1[i].pos);
                    continue;
                }
                a1[i].genotyped = 1; a2[j].genotyped = 1;
                dp1 = qual1->read_data; dp2 = qtmp->read_data;
                for (l = 0, allele_dp = 0; l < n; ++l) {
                    allele_dp += dp1[i + 1];
                    dp1 += qual1->n_alleles;
                }
                if (allele_dp < min_dp) { fprintf(stderr, "Allele %d less than minimum depth\n", i); continue; }
                for (l = 0, allele_dp = 0; l < n; ++l) {
                    allele_dp += dp2[j + 1];
                    dp2 += qtmp->n_alleles;
                }
                if (allele_dp < min_dp) { fprintf(stderr, "Allele %d less than minimum depth\n", i); continue; }
                if (n_var >= m_var) {
                    m_var = m_var ? m_var << 1 : 1;
                    qual2 = (qual_sum_t**)realloc(qual2, m_var * sizeof(qual_sum_t*));
                    var_idx = (int*)realloc(var_idx, m_var * sizeof(int) * 2);
                }
                fprintf(stderr, "At %d:%d, adding i=%d and j=%d\n", qual1->tid, qual1->pos, i, j);
                var_idx[n_var * 2] = i;
                var_idx[n_var * 2 + 1] = j;
                qual2[n_var++] = qtmp;
            }
            if (n_var) { // variants to genotype were found
                double div_sum = 0;
                int mq_sum = 0, qlen_sum = 0, as_sum = 0, n_reads = 0, as_score = 0, as_reads = 0, dp, dp2;
                fprintf(stderr, "%d variants to genotype\n", n_var);
                gts = gts_init(n, n_var);
                fprintf(stderr, "Genotypes initalized\n");
                genotype_samples(gts, qual1, qual2, var_idx, n_var, n);
                fprintf(stderr, "Genotyping finished\n");
                printf("%s\t%d\t.\tN\t", h->target_name[qual1->tid], qual1->pos);
                for (i = 0; i < n_var; ++i) {
                    // print allele information //
                    // do not treat deletions/insertions specially //
                    char dir = qual1->alleles[var_idx[i * 2]].ori2 ? ']' : '[';
                    if (i) putchar(',');
                    if (qual1->alleles[var_idx[i * 2]].ori1) {
                        printf("N%c%s:%d%c", dir, h->target_name[qual2[i]->tid], qual2[i]->pos, dir);
                    } else {
                        printf("%c%s:%d%cN", dir, h->target_name[qual2[i]->tid], qual2[i]->pos, dir);
                    }
                }
                fprintf(stderr, "Alleles printed at %s\t%d\n", h->target_name[qual1->tid], qual1->pos);
                printf("\t.\t.\tQGAP=");
                for (i = 0; i < n_var; ++i) {
                    // print the info field //
                    if (i) putchar(',');
                    printf("%d", qual1->alleles[var_idx[i * 2]].qdist);
                    mq_sum += qual1->mq_sum;
                    mq_sum += qual2[i]->mq_sum;
                    div_sum += qual1->div_sum;
                    div_sum += qual2[i]->div_sum;
                    qlen_sum += qual1->qlen_sum;
                    qlen_sum += qual2[i]->qlen_sum;
                    as_sum += qual1->as_sum;
                    as_sum += qual2[i]->as_sum;
                    n_reads += qual1->n_reads;
                    n_reads += qual2[i]->n_reads;
                }
                fprintf(stderr, "Info printed\n");
                printf(";MAPQ=%.2f;DIV=%.3f;QLEN=%.2f;SC=%.2f;AS_SC=", sqrt((double)mq_sum / n_reads), sqrt(div_sum / n_reads), sqrt((double)qlen_sum / n_reads), sqrt((double)as_sum / n_reads));
                as_score += qual1->as_score[0];
                as_reads += qual1->as_reads[0];
                for (i = 0; i < n_var; ++i) {
                    as_score += qual2[i]->as_score[0];
                    as_reads += qual2[i]->as_reads[0];
                }
                printf("%.0f", sqrt((double)as_score / as_reads));
                as_score = 0; as_reads = 0;
                for (i = 0; i < n_var; ++i) {
                    as_score += qual1->as_score[var_idx[i * 2] + 1];
                    as_score += qual2[i]->as_score[var_idx[i * 2 + 1] + 1];
                    as_reads += qual1->as_reads[var_idx[i * 2] + 1];
                    as_reads += qual2[i]->as_reads[var_idx[i * 2 + 1] + 1];
                    printf(",%.0f", sqrt((double)as_score / as_reads));
                }
                for (dp = 0, i = 0; i < n; ++i) {
                    dp += qual1->read_data[(i * qual1->n_alleles)];
                }
                printf(";BPD=%d", dp);
                for (dp = 0, i = 0; i + 1 < n_var; ++i) {
                    for (j = 0; j < n; ++j) {
                        dp += qual2[i]->read_data[qual2[i]->n_alleles * j];
                    }
                    printf(",%d", dp);
                }
                printf(",%d", dp);
                for (dp = 0, dp2 = 0, i = 0; i < n_var; ++i) {
                    for (j = 0; j < n; ++j) {
                        dp += qual1->read_data[(qual1->n_alleles * j) + var_idx[i * 2] + 1];
                        dp2 += qual2[i]->read_data[(qual2[i]->n_alleles * j) + var_idx[i * 2 + 1] + 1];
                    }
                    printf(",%d,%d", dp, dp2);
                }
                fprintf(stderr, "More info printed\n");
                if (n_var < 2) { // Only analyze HWE and consistency at diploid sites
                    int n_rhom = 0, n_het = 0, n_ahom = 0, res, n_gt;
                    double maf, chi2 = 0, expected, hwe;
                    for (i = 0; i < n; ++i) {
                        if (gts[i].genotype_confidence < 30) continue;
                        if (gts[i].gt == 0) {
                            ++n_rhom;
                        } else if (gts[i].gt == 1) {
                            ++n_het;
                        } else {
                            ++n_ahom;
                        }
                    }
                    n_gt = n_het + n_rhom + n_ahom;
                    maf = ((double)n_ahom * 2.0 + (double)n_het) / (double)((n_rhom + n_het + n_ahom) * 2);
                    fprintf(stderr, "hom ref = %d, het = %d, hom alt = %d, maf = %f\n", n_rhom, n_het, n_ahom, maf);
                    expected = n_gt * (1 - maf) * (1 - maf);
                    if (expected > 0.001) {
                        chi2 += (n_rhom - expected) * (n_rhom - expected) / expected;
                    }
                    expected = n_gt * 2 * maf * (1 - maf);
                    if (expected > 0.001) {
                        chi2 += (n_het - expected) * (n_het - expected) / expected;
                    }
                    expected = n_gt * maf * maf;
                    if (expected > 0.001) {
                        chi2 += (n_ahom - expected) * (n_ahom - expected) / expected;
                    }
                    if (chi2 > 0.0001) {
                        fprintf(stderr, "Finding HWE with %f and %f\n", chi2 / 2, 0.5);
                        hwe = gammds(chi2 / 2, 0.5, &res);
                        fprintf(stderr, "Found HWE=%f with result %d\n", hwe, res);
                    } else {
                        hwe = 0.0000001;
                    }
                    printf(";HWE=%f", log10(1 - hwe) * -10);
                    if (ped_h) {
                        printf(";GT_SEG=%f", log_segregation(gts, n, mi_prob, ped_h));
                    }
                }
                // Genotype //
                fprintf(stderr, "Printing genotype\n");
                print_genotypes(gts, qual1, qual2, var_idx, n_var, n);
                putchar('\n');
                fprintf(stderr, "Destroying genotype array\n");
                gts_destroy(gts);
            }
        }
    }
    free(var_idx);
    free(qual2);
    return 1;
}
