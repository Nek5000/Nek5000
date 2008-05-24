S     = .
L     = 
P     = 

FLAGS = 
OLAGS =
prefix = $(bin_nek_tools)

NOBJ1 = reatore2.o byte.o strings.o
NOBJ2 = re2torea.o byte.o strings.o

all: reatore2

reatore2:	$(NOBJ1)
	$(F77) -o $(prefix)/rea2re2 $(NOBJ1) 

clean:
	'rm' *.o

reatore2.o  	: reatore2.f	; $(F77) -c $P $(FLAGS)  reatore2.f 
re2torea.o  	: re2torea.f	; $(F77) -c $P $(FLAGS)  re2torea.f 
strings.o  	: strings.f	; $(F77) -c $P $(FLAGS)  strings.f 
byte.o		: byte.c	; $(CC)  -c -O2  byte.c
