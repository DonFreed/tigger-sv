#include <stdio.h>
#include "htslib/kstring.h"
#include "ped.h"

khash_t(colmap) *map_samples(char **samples, int n)
{
    int i, ret;
    khash_t(colmap) *h = kh_init(colmap);
    khiter_t k;
    for (i = 0; i < n; ++i) {
        k = kh_put(colmap, h, samples[i], &ret);
        if (ret < 0) {
            fprintf(stderr, "Error during column map creation\n");
            return 0;
        }
        kh_value(h, k) = i;
    }
    return h;
}

khash_t(ped) *read_ped(const char *fnped, khash_t(colmap) *smp_cols)
{
    khash_t(ped) *h = kh_init(ped);
    khiter_t k;
    int res, ln = 0, n, *offsets = 0, max = 0, ret = 0;
    kstring_t line = {0, 0, 0};
    FILE *fp = fopen(fnped, "r");
    while ((res = kgetline(&line, (kgets_func *)fgets, fp)) != EOF) {
        ++ln;
        fprintf(stderr, "Line is %s\n", line.s);
        n = ksplit_core(line.s, '\t', &max, &offsets);
        if (n != 6) {
            fprintf(stderr, "Error: ped file is poorly formatted at line %d. Incorrect number of columns\n", ln);
            return 0;
        }
        if (line.s[offsets[2]] != '0') {
            parents_t tmp;
            int child_col;
            k = kh_get(colmap, smp_cols, line.s + offsets[1]);
            if (k == kh_end(smp_cols)) {
                fprintf(stderr, "Warning. Sample %s is present in ped file but abset from input\n", line.s + offsets[1]);
                line.l = 0;
                continue;
            }
            child_col = kh_value(smp_cols, k);
            k = kh_get(colmap, smp_cols, line.s + offsets[2]);
            if (k == kh_end(smp_cols)) {
                fprintf(stderr, "Warning. Sample %s is present in ped file but abset from input\n", line.s + offsets[2]);
                line.l = 0;
                continue;
            }
            tmp.fa_col = kh_value(smp_cols, k);
            k = kh_get(colmap, smp_cols, line.s + offsets[3]);
            if (k == kh_end(smp_cols)) {
                fprintf(stderr, "Warning. Sample %s is present in ped file but abset from input\n", line.s + offsets[3]);
                line.l = 0;
                continue;
            }
            tmp.mo_col = kh_value(smp_cols, k);
            k = kh_put(ped, h, (khint32_t)child_col, &ret);
            kh_value(h, k) = tmp;
        }
        line.l = 0;
    }
    free(line.s);
    free(offsets);
    return h;
}

