MPI ?= 1
MPIIO ?= 1
ADDUS ?= 1
USREXIT ?= 0
LIBNAME ?= gs
BLAS ?= 0
CFLAGS ?= -O2

ifneq (0,$(MPI))
  G+=-DMPI
endif

ifneq (0,$(MPIIO))
  ifneq (0,$(MPI)) 
    G+=-DUSEMPIIO
  endif
endif

ifneq (0,$(ADDUS))
  G+=-DUNDERSCORE
endif

ifneq (0,$(USREXIT))
  G+=-DUSE_USR_EXIT
endif

ifeq (0,$(BLAS))
  G+=-DUSE_NAIVE_BLAS
endif

ifeq (1,$(BLAS))
  G+=-DUSE_CBLAS
endif

ifneq ($(PREFIX),)
  G+=-DPREFIX=$(PREFIX)
endif

ifneq ($(FPREFIX),)
  G+=-DFPREFIX=$(FPREFIX)
endif

G+=-DGLOBAL_LONG_LONG
#G+=-DPRINT_MALLOCS=1
#G+=-DGS_TIMING -DGS_BARRIER

CCCMD=$(CC) $(CFLAGS) $(G)
LINKCMD=$(CC) $(CFLAGS) $(G) $^ -o $@ -L./ -l$(LIBNAME) -lm $(LDFLAGS)

TESTS=sort_test sort_test2 sarray_sort_test spchol_test \
      comm_test crystal_test sarray_transfer_test \
      gs_test gs_unique_test \
      xxt_test xxt_test2 crs_test \
      findpts_el_2_test findpts_el_2_test2 \
      findpts_el_3_test findpts_el_3_test2 \
      findpts_local_test findpts_test \
      poly_test poly_test2 lob_bnd_test obbox_test

GS=gs.o sort.o sarray_transfer.o sarray_sort.o \
   gs_local.o fail.o crystal.o comm.o tensor.o
XXT=sparse_cholesky.o xxt.o
AMG=amg.o
FWRAPPER=fcrystal.o fcrs.o findpts.o
INTP=findpts_local.o obbox.o poly.o lob_bnd.o findpts_el_3.o findpts_el_2.o

.PHONY: all lib deps tests clean objects odepinfo

all : lib tests

lib: $(GS) $(XXT) $(AMG) $(FWRAPPER) $(INTP) rand_elt_test.o
	@$(AR) cr lib$(LIBNAME).a $?
	@ranlib lib$(LIBNAME).a

tests: $(TESTS)

clean: ; @$(RM) *.o *.s *.a $(TESTS)

cmds: ; @echo CC = $(CCCMD); echo LINK = $(LINKCMD);

deps: ; ./cdep.py *.c > makefile.cdep;

odepinfo: deps objects; @./odep_info.py *.o

$(TESTS): % : %.o | lib
	$(LINKCMD)

-include makefile.cdep

%.o: %.c ; $(CCCMD) -c $< -o $@
%.s: %.c ; $(CCCMD) -S $< -o $@
objects: $(OBJECTS) ;

#poly_imp.h: gen_poly_imp.c
#	$(RM) poly_imp.h;
#	$(CC) -lgmp -lm gen_poly_imp.c -o gen_poly_imp;
#	./gen_poly_imp > poly_imp.h;
#	$(RM) gen_poly_imp
