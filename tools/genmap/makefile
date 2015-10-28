FLAGS =  $(FFLAGS_IN)
CFLAGS =  $(CFLAGS_IN)
LFLAGS = $(LFLAGS_IN)
prefix = $(bin_nek_tools)


OBJS = genmap.o byte.o 

all: genmap

genmap:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/genmap $(OBJS) $(LFLAGS)

clean:
	'rm' *.o

genmap.o          : genmap.f		 		;  $(F77) -c $(FLAGS) genmap.f
byte.o				: ../../nek/byte.c	;  $(CC)  -c $(CFLAGS) ../../nek/byte.c
