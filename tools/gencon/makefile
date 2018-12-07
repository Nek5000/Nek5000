prefix = $(bin_nek_tools)

OBJS = gencon.o byte.o 

all: gencon

gencon:	$(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/gencon $(OBJS) $(LDFLAGS)

clean:
	@rm -f *.o

gencon.o        : gencon.f		;  $(FC) -c $(FFLAGS) gencon.f
byte.o		: ../../core/byte.c	;  $(CC) -c $(CFLAGS) ../../core/byte.c
