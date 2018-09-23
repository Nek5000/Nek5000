prefix = $(bin_nek_tools)

OBJS = genmap.o byte.o 

all: genmap

genmap:	$(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/genmap $(OBJS) $(LDFLAGS)

clean:
	@rm -f *.o

genmap.o        : genmap.f		;  $(FC) -c $(FFLAGS) genmap.f
byte.o		: ../../core/byte.c	;  $(CC) -c $(CFLAGS) ../../core/byte.c
