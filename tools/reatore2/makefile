FFLAGS = 
CFLAGS = -DUNDERSCORE
prefix = $(bin_nek_tools)

NOBJ1 = reatore2.o byte.o strings.o
NOBJ2 = re2torea.o byte.o strings.o

all: reatore2

reatore2:	$(NOBJ1)
	$(F77) -o $(prefix)/reatore2 $(NOBJ1) 

clean:
	'rm' *.o

reatore2.o  	: reatore2.f		 ; $(F77) -c $(FFLAGS)  reatore2.f 
re2torea.o  	: re2torea.f	 	 ; $(F77) -c $(FFLAGS)  re2torea.f 
strings.o  	   : strings.f			 ; $(F77) -c $(FFLAGS)  strings.f 
byte.o		   : ../../nek/byte.c  ; $(CC)  -c $(CFLAGS) ../../nek/byte.c
