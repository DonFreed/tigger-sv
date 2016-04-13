#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include "tigger_version.h"
#include "genotype.h"

void print_header(bam_hdr_t *h, int optind, int n, char *argv[])
{
    int i;
    printf("##fileformat=VCFv4.2\n");
    printf("##source=tigger-sv-%s\n", TIGGERSV_VERSION);
    printf("##ALT=<ID=DEL,Description=\"Deletion relative to the reference genome\">\n");
    printf("##ALT=<ID=INS,Description=\"Insertion relative to the reference genome\">\n");
    printf("##ALT=<ID=INV,Description=\"Inversion relative to the reference genome\">\n");
    printf("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the structural variant along the reference\">\n");
    printf("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
    printf("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of the variant\">\n");
    printf("##INFO=<ID=QGAP,Number=1,Type=Integer,Description=\"Length of the gap on the query sequence\">\n");
    printf("##INFO=<ID=MAPQ,Number=1,Type=Float,Description=\"RMS mapping quality of supporing alignments\">\n");
    printf("##INFO=<ID=DIV,Number=1,Type=Float,Description=\"RMS divergence of supporting alignments\">\n");
    //printf("##INFO=<ID=SC,Number=1,Type=Integer,Description=\"RMS alignment score of supporting alignments\">\n");
    //printf("##INFO=<ID=TIPQ,Number=1,Type=Integer,Description=\"RMS quality/depth of the supporting alignments\">\n");
    for (i = 0; i < h->n_targets; ++i) {
        printf("##contig=<ID=%s,length=%" PRIu32 ">\n", h->target_name[i], h->target_len[i]);
    }
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (i = 0; i < n; ++i) {
        printf("\t%s", argv[optind + i]);
    }
    putchar('\n');
}

int genotype_sv(bam_hdr_t *h, int n, khash_t(sv_geno) *geno_h)
{
    int i, j;
    khiter_t k1, k2;
    qual_sum_t qual1, qual2;
    allele_t *a1, *a2;
    uint16_t *dp;
    uint64_t id;
    for (k1 = kh_begin(geno_h); k1 != kh_end(geno_h); ++k1) {
        if (kh_exist(geno_h, k1)) {
            qual1 = kh_value(geno_h, k1);
            a1 = qual1.alleles;
            for (i = 0; i < qual1.n_alleles; ++i) {
                int tid, pos, end_tid, end_pos, ori1, ori2;
                if (a1[i].genotyped) { continue; }
                id = (uint64_t)a1[i].tid << 32 | a1[i].pos;
                k2 = kh_get(sv_geno, geno_h, id);
                if (k2 == kh_end(geno_h)) {
                    fprintf(stderr, "Error: breakpoint mate not found.\n");
                    return -1;
                }
                qual2 = kh_value(geno_h, k2);
                a2 = qual2.alleles;
                for (j = 0; j < qual2.n_alleles; ++j) {
                    if (qual1.tid == a2[j].tid && qual1.pos == a2[j].pos) break;
                }
                if (j == qual2.n_alleles) {
                    fprintf(stderr, "Error: breakpoint does not have a matching allele.\n");
                    return -1;
                }
                a1[i].genotyped = 1; a2[j].genotyped = 1;
                if (qual1.tid < qual2.tid || (qual1.tid == qual2.tid && qual1.pos < qual2.pos)) {
                    tid = qual1.tid;
                    pos = qual1.pos;
                    end_tid = qual2.tid;
                    end_pos = qual2.pos;
                    ori1 = a1[i].ori1;
                    ori2 = a1[i].ori2;
                } else {
                    tid = qual2.tid;
                    pos = qual2.pos;
                    end_tid = qual1.tid;
                    end_pos = qual1.pos;
                    ori1 = a1[i].ori2;
                    ori1 = a1[i].ori1;
                }
                printf("%s\t%d\t.\tN", h->target_name[tid], pos);
                if (a1[i].type == 'D') {
                    printf("\t<DEL>\t.\tSVTYPE=DEL;END=%d;", end_pos);
                } else if (a1[i].type == 'I') {
                    printf("\t<INS>\t.\tSVTYPE=INS;END=%d;", end_pos);
                } else {
                    char dir = ori2 ? ']' : '[';
                    if (ori1) {
                        printf("\tN%c%s:%d%c\t.\t", dir, h->target_name[end_tid], end_pos, dir);
                    } else {
                        printf("\t%c%s:%d%cN\t.\t", dir, h->target_name[end_tid], end_pos, dir);
                    }
                }
                printf("QGAP=%d;", a1[i].qdist);
                printf("MAPQ=%.2f;", sqrt((double)(qual1.mq_sum + qual2.mq_sum) / (qual1.n_reads + qual2.n_reads)));
                printf("DIV=%.3f;", sqrt((qual1.div_sum + qual2.div_sum) / (qual1.n_reads + qual2.n_reads)));
                printf("QLEN=%.2f\n", sqrt((double)(qual1.qlen_sum + qual2.qlen_sum) / (qual1.n_reads + qual2.n_reads)));
                
            }
        }
    }
    return 1;
}
