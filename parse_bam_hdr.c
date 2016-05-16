#include <string.h>
#include <stdlib.h>
#include "htslib/kstring.h"
#include "parse_bam_hdr.h"

char *find_sample(bam_hdr_t *hdr, int *res)
{
    kstring_t s = {0, 0, 0};
    int max = 0, *offsets = 0, i, n;
    int line_max = 0, *line_off = 0, line_n, j, str_len;
    char *sample = 0;
    kputsn(hdr->text, (int)hdr->l_text, &s);
    n = ksplit_core(s.s, '\n', &max, &offsets);
    *res = 0;
    for (i = 0; i < n; ++i) {
        if (s.s[offsets[i]] == '@' && s.s[offsets[i] + 1] == 'R' && s.s[offsets[i] + 2] == 'G') {
            line_n = ksplit_core(s.s + offsets[i], '\t', &line_max, &line_off);
            for (j = 0; j < line_n; ++j) {
                if (s.s[offsets[i] + line_off[j]] == 'S' && s.s[offsets[i] + line_off[j] + 1] == 'M' && s.s[offsets[i] + line_off[j] + 2] == ':') {
                    str_len = strlen(s.s + offsets[i] + line_off[j] + 3);
                    if (sample) {
                        if ((strcmp(sample, s.s + offsets[i] + line_off[j] + 3)) != 0) {
                            fprintf(stderr, "Error. Multiple samples %s and %s found in a single bam file.\n", sample, s.s + offsets[i] + line_off[j] + 3);
                            *res = -1;
                        }
                    } else {
                        sample = (char*)malloc(str_len + 1);
                        strcpy(sample, s.s + offsets[i] + line_off[j] + 3);
                        sample[str_len] = '\0';
                    }
                }
            }
        } 
    }
    
    free(s.s);
    free(offsets);
    free(line_off);
    if (*res) {
        return 0;
    } else {
        return sample;
    }
}

