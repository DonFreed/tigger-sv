#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "plp2sv.h"

typedef struct {
    samFile *fp;
    bam_hdr_t *hdr;
    int min_mapq, min_as, min_len;
} aux_t;

int bam_cigar2ulen(int n_cigar, const uint32_t *cigar)
{
    int k, l;
    for (k = l = 0; k < n_cigar; ++k) {
        if (bam_cigar_type(bam_cigar_op(cigar[k])) &1) {
            l += bam_cigar_oplen(cigar[k]);
        } else if (bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP) {
            l += bam_cigar_oplen(cigar[k]);
        }
    }
    return l;
}

int read_bam(void *data, bam1_t *b)
{
    aux_t *aux = (aux_t*)data;
    int ret;
    while(1) {
        uint8_t *tmp = 0;
        ret = sam_read1(aux->fp, aux->hdr, b);
        if (ret < 0) break;
        if (b->core.flag & (BAM_FUNMAP)) continue;
        if ((int)b->core.qual < aux->min_mapq) continue;
        if (bam_cigar2ulen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len) continue;
        tmp = bam_aux_get(b, "AS");
        if (tmp && bam_aux2i(tmp) < aux->min_as) continue;
        break;
    }
    return ret;
}

typedef struct {
    int min_q, min_s, min_len, min_dp;
} cmdopt_t;

void usage(FILE *fp, cmdopt_t *o)
{
    
    fprintf(fp, "\n");
    fprintf(fp, "tigger-sv\n");
    fprintf(fp, "  Identify structual variants from unitigs in cohorts of samples\n");
    fprintf(fp, "\n");
    fprintf(fp, "Usage: tigger-sv [options] <in1.bam> [in2.bam ...]\n");
    fprintf(fp, "\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -h           print this message.\n");
    fprintf(fp, "  -q INT       minimum mapping quality for alignments [%d]\n", o->min_q);
    fprintf(fp, "  -s INT       minimum alignment score [%d]\n", o->min_s);
    fprintf(fp, "  -l INT       minimum unitig length [%d]\n", o->min_len);
    fprintf(fp, "  -d INT       minimum breakpoint depth [%d]\n\n", o->min_dp);
}

int main(int argc, char *argv[])
{
    int c, i, n, ret;
    int tid, pos, *n_plp;
    cmdopt_t o;
    bam_mplp_t mplp;
    const bam_pileup1_t **plp;
    aux_t **data;
    bam_hdr_t *h = 0;
    sv_vec_t sv = {0, 0, 0};
    sv_t sv1;
    khiter_t k_iter;
    khash_t(sv_hash) *sv_h = kh_init(sv_hash);
    
    o.min_q = 40; o.min_s = 80; o.min_len = 150; o.min_dp = 10;
    while ((c = getopt(argc, argv, "hq:s:l:d:")) >= 0) {
        if (c == 'h') { usage(stderr, &o); return 0; }
        else if (c == 'q') o.min_q = atoi(optarg);
        else if (c == 's') o.min_s = atoi(optarg);
        else if (c == 'l') o.min_len = atoi(optarg);
        else if (c == 'd') o.min_dp = atoi(optarg);
    }
    
    if (argc - optind < 1) {
        usage(stderr, &o);
        return 1;
    }

    // Open files and initalize aux data //
    n = argc - optind;
    data = calloc(n, sizeof(aux_t*));
    for (i = 0; i < n; ++i) {
        data[i] = calloc(1, sizeof (aux_t));
        data[i]->fp = sam_open(argv[optind + i], "r");
        if (!data[i]->fp) {
            fprintf(stderr, "Input file \"%s\" could not be opened.\n", argv[optind + 1]);
            return 1;
        }
        data[i]->min_mapq = o.min_q;
        data[i]->min_as = o.min_s;
        data[i]->min_len = o.min_len;
        data[i]->hdr = sam_hdr_read(data[i]->fp);
        if (!data[i]->hdr) {
            fprintf(stderr, "Could not read the header for input file \"%s\".\n", argv[optind + 1]);
            return 1;
        }
    }
    h = data[0]->hdr;

    // The core data processing loop //
    mplp = bam_mplp_init(n, read_bam, (void**)data);
    n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads in mplp
    while ((ret = bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // iterate of positions with coverage
        int n_sv;
        sv.n = 0;
        n_sv = plp2sv(h, tid, pos, n, n_plp, plp, sv_h);
        if (n_sv) { fprintf(stderr, "SV detected at %d:%d\n", tid, pos); }
    }
    for (k_iter = kh_begin(sv_h); k_iter != kh_end(sv_h); ++k_iter) {
        if (kh_exist(sv_h, k_iter)) {
            sv1 = kh_value(sv_h, k_iter);
            fprintf(stderr, "SV tid1=%d, tid2=%d, pos1=%d, pos2=%d, ori1=%d, ori2=%d\n", sv1.tid1, sv1.tid2, sv1.pos1, sv1.pos2, sv1.ori1, sv1.ori2);
        }
    }


    free(n_plp);
    free(plp);
    bam_mplp_destroy(mplp);
    for (i = 0; i < n; ++i) { 
        bam_hdr_destroy(data[i]->hdr);
        sam_close(data[i]->fp);
        free(data[i]); 
    }
    free(data);
    free(sv.sv);
    kh_destroy(sv_hash, sv_h);
    return 0;
}
