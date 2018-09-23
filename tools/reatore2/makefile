prefix = $(bin_nek_tools)

NOBJ1 = reatore2.o byte.o strings.o
NOBJ2 = re2torea.o byte.o strings.o

all: reatore2 re2torea

reatore2:	$(NOBJ1)
	$(FC) -o $(prefix)/reatore2 $(NOBJ1) 

re2torea:	$(NOBJ2)
	$(FC) -o $(prefix)/re2torea $(NOBJ2) 


clean:
	@rm -f *.o

reatore2.o  	: reatore2.f		 ; $(FC) -c $(FFLAGS)  reatore2.f 
re2torea.o  	: re2torea.f	 	 ; $(FC) -c $(FFLAGS)  re2torea.f 
strings.o  	: strings.f		 ; $(FC) -c $(FFLAGS)  strings.f 
byte.o		: byte.c     		 ; $(CC)  -c $(CFLAGS) byte.c
