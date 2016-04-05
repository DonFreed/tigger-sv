CC=				gcc
CFLAGS=			-g -Wall
CPPFLAGS=		-Ihtslib
OBJS=			tigger-sv.o plp2sv.o cigar.o
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


tigger-sv.o: plp2sv.h array.h
plp2sv.o: cigar.h array.h
