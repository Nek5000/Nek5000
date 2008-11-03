S     = .
L     = 
P     = -r8
F77   = pgf77
CC    = pgcc

FLAGS =  -O2  -mcmodel=medium # -g
OLAGS =

NOBJ1 = nekmerge.o reader.o byte.o strings.o tsort.o

all: nekmerge

nekmerge:	$(NOBJ1)
	$(F77) -o nekmerge $(NOBJ1) $(FLAGS)

clean:
	'rm' *.o

nekmerge.o  	: nekmerge.f	; $(F77) -c $P $(FLAGS)  nekmerge.f 
reader.o  	: reader.f	; $(F77) -c $P $(FLAGS)  reader.f 
tsort.o  	: tsort.f	; $(F77) -c $P $(FLAGS)  tsort.f 
strings.o  	: strings.f	; $(F77) -c $P $(FLAGS)  strings.f 
byte.o		: byte.c	; $(CC)  -c -O2  byte.c
