FLAGS =  $(FFLAGS_IN)
LFLAGS = $(LFLAGS_IN)
prefix = $(bin_nek_tools)

OBJS = genbox.o 

all: genbox


genBox:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/genbox $(OBJS) $(LFLAGS) 

clean:
	'rm' *.o

genbox.o          : genbox.f		 ;  $(F77) -c $(FLAGS) genbox.f 
