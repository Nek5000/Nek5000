FLAGS = 
LFLAGS = 
prefix = $(bin_nek_tools)

OBJS = genBox.o 

all: genBox


genBox:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/genBox $(OBJS) $(LFLAGS) 

clean:
	'rm' *.o

genBox.o          : genBox.f		 ;  $(F77) -c $(FLAGS) genBox.f 
