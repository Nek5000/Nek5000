FFLAGS = 
CFLAGS = -DUNDERSCORE 

prefix = $(bin_nek_tools)

OBJS = bin2asc.o byte.o subs.o

all: bin2asc


bin2asc:	$(OBJS)
	$(F77) $(FFLAGS) -o $(prefix)/bin2asc $(OBJS) 

clean:
	'rm' *.o 


bin2asc.o        : bin2asc.f	   	       ;  $(F77) -c $(FFLAGS) bin2asc.f
subs.o           : subs.f 					  	 ;  $(F77) -c $(FFLAGS) subs.f
byte.o           : ../../nek/byte.c        ;  $(CC) -c $(CFLAGS) ../../nek/byte.c
