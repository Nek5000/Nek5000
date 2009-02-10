FLAGS =  $(FFLAGS_IN)
LFLAGS = 
prefix = $(bin_nek_tools)

OBJS = genbox.o 

all: genbox


genbox:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/genbox $(OBJS) $(LFLAGS) 

clean:
	'rm' *.o

genbox.o          : genbox.f		 ;  $(F77) -c $(FLAGS) genbox.f 
