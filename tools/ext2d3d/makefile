FLAGS = 
LFLAGS = 
prefix = $(bin_nek_tools)

OBJS = n2to3.o 

all: ext2d3d


ext2d3d:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/ext2d3d $(OBJS) $(LFLAGS) 

clean:
	'rm' *.o

n2to3.o          : n2to3.f		 ;  $(F77) -c $(FLAGS) n2to3.f 
