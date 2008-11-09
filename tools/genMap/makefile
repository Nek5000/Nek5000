FLAGS = -r8 
LFLAGS = 
prefix = $(bin_nek_tools)

OBJS = genMap.o 

all: genMap


genMap:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/genMap $(OBJS) $(LFLAGS) 

clean:
	'rm' *.o

genMap.o          : genMap.f		 ;  $(F77) -c $(FLAGS) genMap.f 
