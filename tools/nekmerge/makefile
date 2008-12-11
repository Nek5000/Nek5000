S     = .
L     = 
P     = -r8

FLAGS =  -DUNDERSCORE # -g
OLAGS =
prefix = $(bin_nek_tools)

NOBJ1 = nekmerge.o reader.o byte.o strings.o tsort.o

all: nekmerge

nekmerge:	$(NOBJ1)
	$(F77) -o  $(prefix)/nekmerge $(NOBJ1) $(FLAGS)

clean:
	'rm' *.o

nekmerge.o	: nekmerge.f; $(F77) -c $P $(FLAGS)  nekmerge.f 
reader.o  	: reader.f	; $(F77) -c $P $(FLAGS)  reader.f 
tsort.o  	: tsort.f	; $(F77) -c $P $(FLAGS)  tsort.f 
strings.o  	: strings.f	; $(F77) -c $P $(FLAGS)  strings.f 
byte.o		: ../../nek/byte.c	; $(CC) -c $(FLAGS)  ../../nek/byte.c
