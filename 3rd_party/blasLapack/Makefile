AR = ar
SDIR = ./
ARFLAGS = -cr

.SUFFIXES: .f

SRC = $(wildcard $(SDIR)/*f)
OBJS = $(subst $(SDIR)/,$(SDIR),$(SRC:.f=.o))

.f.o:
	$(FC) $(FFLAGS) -c $*.f 

lib: $(OBJS)
	$(AR) $(ARFLAGS) libblasLapack.a *.o

clean:
	@rm -f libblasLapack.a *.o
