FLAGS =  $(FFLAGS_IN)
CFLAGS = $(CFLAGS_IN)
LFLAGS = $(LFLAGS_IN) -lexodus -lnetcdf -lexodus_for
prefix = $(bin_nek_tools)


OBJS = exo2nek.o byte.o speclib.o mxm.o 

all: exo2nek

exo2nek:$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/exo2nek $(OBJS) $(LFLAGS)

clean:
	'rm' *.o

exo2nek.o	: exo2nek.f			;  $(F77) -c $(FLAGS) exo2nek.f
byte.o		: ../../core/byte.c		;  $(CC)  -c $(CFLAGS) ../../core/byte.c
mxm_std.o	: ../../core/mxm_std.f		;  $(F77) -c $(FLAGS) ../../core/mxm_std.f
speclib.o	: ../../core/speclib.f		;  $(F77) -c $(FLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(F77) -c $(FLAGS) mxm.f
