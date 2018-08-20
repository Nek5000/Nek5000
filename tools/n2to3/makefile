prefix = $(bin_nek_tools)

OBJS = n2to3.o byte.o

all: n2to3


n2to3:	$(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/n2to3 $(OBJS) $(LDFLAGS) 

clean:
	@rm -f *.o

n2to3.o          : n2to3.f		 ;  $(FC) -c $(FFLAGS) $(R8) n2to3.f 
byte.o           : byte.c      		 ;  $(CC)  -c $(CFLAGS) byte.c
