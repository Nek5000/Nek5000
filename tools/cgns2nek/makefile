prefix = $(bin_nek_tools)
INC = -I ./cgns/include
LIBS = -L ./cgns/lib -lcgns

OBJS = module_global.o write_re2_file.o read_cgns.o cgns2nek.o byte.o  

all: lib cgns2nek

cgns2nek: $(OBJS)
	 $(FC) $(FFLAGS) -o $(prefix)/cgns2nek $^ $(LIBS) $(LDFLAGS)

clean:
	@rm -f *.o *.mod 
	@cd ./cgns ; ./install clean

lib:
	@cd cgns; env CC="$(CC)" CFLAGS="$(BIGMEM)" FC="$(FC)" FFLAGS="$(BIGMEM)" ./install

t.o  		 : t.F90			;  $(FC) -c $(FFLAGS) $(INC) t.F90 
module_global.o  : module_global.f90		;  $(FC) -c $(FFLAGS) module_global.f90 
write_re2_file.o : write_re2_file.f90		;  $(FC) -c $(FFLAGS) $(INC) write_re2_file.f90
read_cgns.o      : read_cgns.f90		;  $(FC) -c $(FFLAGS) $(INC) read_cgns.f90
cgns2nek.o       : cgns2nek.f90     		;  $(FC) -c $(FFLAGS) $(INC) cgns2nek.f90
byte.o		 : ../../core/byte.c		;  $(CC) -c $(CFLAGS) ../../core/byte.c
