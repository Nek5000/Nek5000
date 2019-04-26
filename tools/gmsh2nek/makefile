prefix = $(bin_nek_tools)

OBJS = mod_SIZE.o gmsh2nek.o byte.o speclib.o mxm.o 

all: gmsh2nek

gmsh2nek: $(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/gmsh2nek $^ $(LDFLAGS)

clean:
	@rm -f *.o 

mod_SIZE.o  : mod_SIZE.f90          ;  $(FC) -c $(FFLAGS) mod_SIZE.f90
gmsh2nek.o	: gmsh2nek.f90			;  $(FC) -c $(FFLAGS) gmsh2nek.f90
byte.o		: ../../core/byte.c		;  $(CC) -c $(CFLAGS) ../../core/byte.c
speclib.o	: ../../core/speclib.f		;  $(FC) -c $(FFLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(FC) -c $(FFLAGS) mxm.f
