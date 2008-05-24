FLAGS =  
prefix = $(bin_nek_tools)

OBJS = bin2asc.o byte.o subs.o

all: bin2asc


bin2asc:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/bin2asc $(OBJS) 

clean:
	'rm' *.o 


bin2asc.o        : bin2asc.f	   	       ;  $(F77) -c $(FLAGS) bin2asc.f
subs.o           : ../interpolation/subs.f ;  $(F77) -c $(FLAGS) ../interpolation/subs.f
byte.o           : ../nek2tec/byte.c       ;  $(CC) -c ../nek2tec/byte.c
