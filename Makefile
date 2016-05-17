CC=				gcc
CFLAGS=			-g -Wall -O2
CPPFLAGS=		-Ihtslib
OBJS=			tigger-sv.o plp2sv.o cigar.o sv_qual.o mempool.o genotype.o bedidx.o asa147.o parse_bam_hdr.o ped.o
PROG=			tigger-sv

.SUFFIXES:.c .o
.PHONY:all lib

.c.o:
				$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

all:$(PROG)

lib:
				cd htslib; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhts.a || exit 1; cd ..

tigger-sv:lib $(OBJS)
				$(CC) $(CFLAGS) -o $@ $(OBJS) htslib/libhts.a -lpthread -lz -lm

clean:
				rm -rf gmon.out *.o a.out *~ $(PROG); cd htslib; $(MAKE) clean; cd ..


tigger-sv.o: plp2sv.h array.h mempool.h sv_qual.h genotype.h parse_bam_hdr.h ped.h
plp2sv.o: cigar.h array.h mempool.h
sv_qual.o: cigar.h plp2sv.h array.h
genotype.o: sv_qual.h asa147.h
