prefix = $(bin_nek_tools)
INC = -I./3rd_party/seacas-exodus/packages/seacas/libraries/exodus_for/include
LIBS += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exoIIv2for32/libexoIIv2for32.a
LIBS += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exodus/libexodus.a
LIBS += ./3rd_party/netcdf/install/lib/libnetcdf.a

OBJS = mod_SIZE.o exo2nek.o byte.o speclib.o mxm.o 

all: lib exo2nek

exo2nek: $(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/exo2nek $^ $(LIBS) $(LDFLAGS)

clean:
	@rm -f *.o *.mod 
	@cd ./3rd_party ; rm -rf seacas-exodus netcdf

lib:
	@cd 3rd_party; env CC="$(CC)" CFLAGS="$(BIGMEM)" FC="$(FC)" FFLAGS="$(BIGMEM)" ./install

mod_SIZE.o      : mod_SIZE.f90			;  $(FC) -c $(FFLAGS) mod_SIZE.f90
exo2nek.o	: exo2nek.f90			;  $(FC) -c $(FFLAGS) $(INC) exo2nek.f90
byte.o		: ../../core/byte.c		;  $(CC) -c $(CFLAGS) ../../core/byte.c
speclib.o	: ../../core/speclib.f		;  $(FC) -c $(FFLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(FC) -c $(FFLAGS) mxm.f
