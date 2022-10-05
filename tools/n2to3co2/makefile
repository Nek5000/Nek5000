prefix = $(bin_nek_tools)

OBJS = n2to3co2.o byte.o

all: n2to3co2


n2to3co2:	$(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/n2to3co2 $(OBJS) $(LDFLAGS) 

clean:
	@rm -f *.o

n2to3co2.o      : n2to3co2.f		;  $(FC) -c $(FFLAGS) n2to3co2.f 
byte.o		: ../../core/byte.c	;  $(CC) -c $(CFLAGS) ../../core/byte.c
