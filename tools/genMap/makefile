FLAGS = 
LFLAGS = 
prefix = $(bin_nek_tools)

OBJS = genMapMetis.o 

all: genMap


genMap:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/genMap $(OBJS) $(LFLAGS) 

clean:
	'rm' *.o

genMapMetis.o          : genMapMetis.f		 ;  $(F77) -c $(FLAGS) genMapMetis.f 
