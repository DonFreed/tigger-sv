#ifndef TGGHDR
#define TGGHDR

#include "htslib/sam.h"

char *find_sample(bam_hdr_t *hdr, int *res);

#endif
