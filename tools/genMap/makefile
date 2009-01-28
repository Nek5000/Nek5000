FFLAGS = -r8 -Ktrap=fp -mcmodel=medium -Mprof=func # -g -O0
CFLAGS = -DUNDERSCORE 

F77 = pgf77
cc  = pgcc 

OBJS = genmap.o byte.o 

all: genmap


genmap:	$(OBJS)
	$(F77) $(FFLAGS) -o genmap $(OBJS)

clean:
	'rm' *.o

genmap.o          : genmap.f		 		;  $(F77) -c $(FFLAGS) genmap.f
byte.o				: ../../nek/byte.c	;  $(CC)  -c $(CFLAGS) ../../nek/byte.c
