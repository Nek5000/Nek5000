prefix = $(bin_nek_tools)

NOBJ1 = nekmerge.o reader.o byte.o strings.o tsort.o

all: nekmerge

nekmerge:	$(NOBJ1)
	$(FC) $(FFLAGS) -o  $(prefix)/nekmerge $(NOBJ1) $(LDFLAGS)

clean:
	@rm -f *.o

nekmerge.o	: nekmerge.f    ; $(FC) -c $(FFLAGS)  nekmerge.f 
reader.o  	: reader.f	; $(FC) -c $(FFLAGS)  reader.f 
tsort.o  	: tsort.f	; $(FC) -c $(FFLAGS)  tsort.f 
strings.o  	: strings.f	; $(FC) -c $(FFLAGS)  strings.f 
byte.o		: byte.c	; $(CC) -c $(CFLAGS)  byte.c
