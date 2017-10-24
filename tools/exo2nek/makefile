prefix = $(bin_nek_tools)
INC = -I./3rd_party/seacas-exodus/packages/seacas/libraries/exodus_for/include
LIBS += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exoIIv2for32/libexoIIv2for32.a
LIBS += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exodus/libexodus.a
LIBS += ./3rd_party/netcdf/install/lib/libnetcdf.a

OBJS = exo2nek.o byte.o speclib.o mxm.o 

all: lib exo2nek

exo2nek: $(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/exo2nek $^ $(LIBS) $(LDFLAGS)

clean:
	@rm -f *.o 
	@cd ./3rd_party ; rm -rf seacas-exo* netcdf* *.tar* *.zip*

lib:
	@cd 3rd_party; env CC="$(CC)" CFLAGS="$(BIGMEM)" FC="$(FC)" FFLAGS="$(BIGMEM)" ./install

exo2nek.o	: exo2nek.f			;  $(FC) -c $(FFLAGS) $(INC) exo2nek.f
byte.o		: ../../core/byte.c		;  $(CC) -c $(CFLAGS) ../../core/byte.c
speclib.o	: ../../core/speclib.f		;  $(FC) -c $(FFLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(FC) -c $(FFLAGS) mxm.f
