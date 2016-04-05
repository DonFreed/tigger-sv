#include <stdio.h>
#include <stdlib.h>

/* Converts a string representation of a CIGAR to sam/bam format */
/*
 * 
 * 
 */ 
int str2cigar(const char *s, int *max, uint32_t *cigar)
{
    
