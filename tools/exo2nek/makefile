FLAGS =  $(FFLAGS_IN) -I./3rd_party/seacas-exodus/packages/seacas/libraries/exodus_for/include
CFLAGS = $(CFLAGS_IN)
LFLAGS = $(LFLAGS_IN)
prefix = $(bin_nek_tools)
LIBS   = ./3rd_party/netcdf/lib/libnetcdf.a
LIBS  += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exodus/libexodus.a 
LIBS  += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exoIIv2for32/libexoIIv2for32.a 

OBJS = exo2nek.o byte.o speclib.o mxm.o 

all: exo2nek

exo2nek: $(LIBS) $(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/exo2nek $^ $(LFLAGS)

clean:
	@rm -f *.o 
	@cd ./3rd_party ; rm -rf seacas-exo* netcdf* *.tar* *.zip*

%.a:
	@cd 3rd_party; ./install

exo2nek.o	: exo2nek.f			;  $(F77) -c $(FLAGS) exo2nek.f
byte.o		: ../../core/byte.c		;  $(CC)  -c $(CFLAGS) ../../core/byte.c
speclib.o	: ../../core/speclib.f		;  $(F77) -c $(FLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(F77) -c $(FLAGS) mxm.f
