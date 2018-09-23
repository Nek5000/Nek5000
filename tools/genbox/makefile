prefix = $(bin_nek_tools)

OBJS = genbox.o byte.o 

all: genbox


genbox:	$(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/genbox $(OBJS) $(LDFLAGS) 

clean:
	@rm -f *.o

genbox.o       : genbox.f		;  $(FC) -c $(FFLAGS) genbox.f
byte.o         : byte.c  		;  $(CC)  -c $(CFLAGS) byte.c
