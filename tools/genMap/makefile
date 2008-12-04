FFLAGS = -r8 
CFLAGS = -DUNDERSCORE 
prefix = $(bin_nek_tools)

OBJS = genmap.o byte.o 

all: genmap


genmap:	$(OBJS)
	$(F77) $(FFLAGS) -o $(prefix)/genmap $(OBJS)

clean:
	'rm' *.o

genmap.o          : genmap.f		 		;  $(F77) -c $(FFLAGS) genmap.f
byte.o				: ../../nek/byte.c	;  $(CC)  -c $(CFLAGS) ../../nek/byte.c
