MPI ?= 0
ADDUS ?= 0
LIBNAME ?= gsserial
PREFIX ?= gslib
CBLAS ?= 0

ifneq (0,$(MPI))
  LIBNAME = gs
  G+=-DMPI
endif

ifneq (0,$(ADDUS))
  G+=-DUNDERSCORE
endif

ifneq (0,$(CBLAS))
  G+=-DUSE_CBLAS
endif

G+=-DPREFIX=$(PREFIX)

G+=-DGLOBAL_LONG_LONG
G+=-DAMG_BLOCK_ROWS=2400

G+=-DNEK_EXITT

#G+=-DPRINT_MALLOCS=1
#G+=-DGS_TIMING -DGS_BARRIER
#G+=-DUSE_NAIVE_BLAS

CCCMD=$(CC) -O2 $(G) $(CFLAGS)
LINKCMD=$(CC) $(G) -L./ -l$(LIBNAME) $(LDFLAGS)

TESTS=sort_test sort_test2 sarray_sort_test spchol_test \
      comm_test crystal_test sarray_transfer_test \
      gs_test gs_test_old gs_unique_test \
      xxt_test xxt_test2 crs_test \
#      findpts_el_2_test findpts_el_2_test2 \
#      findpts_el_3_test findpts_el_3_test2 \
#      findpts_local_test findpts_test
#      poly_test poly_test2 lob_bnd_test obbox_test \

GS=gs.o sort.o sarray_transfer.o sarray_sort.o \
   gs_local.o fail.o crystal.o comm.o tensor.o
XXT=sparse_cholesky.o xxt.o
AMG=amg.o
FWRAPPER=fcrystal.o fcrs.o findpts.o
INTP=findpts_local.o obbox.o poly.o lob_bnd.o findpts_el_3.o findpts_el_2.o

.PHONY: all lib cmds deps tests clean objects odepinfo

all : lib tests

lib: $(GS) $(XXT) $(AMG) $(FWRAPPER) $(INTP)
	ar cr lib$(LIBNAME).a $?
	ranlib lib$(LIBNAME).a

tests: $(TESTS) ;

clean: ; $(RM) *.o *.s *.a $(TESTS)

cmds: ; @echo CC = $(CCCMD); echo LINK = $(LINKCMD);

deps: ; ./cdep.py *.c > makefile.cdep;

odepinfo: deps objects; @./odep_info.py *.o

-include makefile.cdep

%.o: %.c ; $(CCCMD) -c $< -o $@
%.s: %.c ; $(CCCMD) -S $< -o $@
objects: $(OBJECTS) ;

#poly_imp.h: gen_poly_imp.c
#	$(RM) poly_imp.h;
#	$(CC) -lgmp -lm gen_poly_imp.c -o gen_poly_imp;
#	./gen_poly_imp > poly_imp.h;
#	$(RM) gen_poly_imp

CRS=$(AMG)

sort_test: sort_test.o ; $(LINKCMD) $^ -o $@
sort_test2: sort_test2.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
sarray_sort_test: sarray_sort_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
spchol_test: spchol_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
comm_test: comm_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
crystal_test: crystal_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
sarray_transfer_test: sarray_transfer_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@

gs_test: gs_test.o $(GS_OBJECTS);		@echo LINK $@; $(LINKCMD) $^ -o $@
gs_test_old: gs_test_old.o $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
gs_unique_test: gs_unique_test.o $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
xxt_test: xxt_test.o $(CRS) $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
xxt_test2: xxt_test2.o $(CRS) $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
crs_test: crs_test.o $(CRS) $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@

poly_test2: poly_test2.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
poly_test: poly_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
lob_bnd_test: lob_bnd_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
obbox_test: rand_elt_test.o obbox_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
findpts_el_2_test2: rand_elt_test.o findpts_el_2_test2.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
findpts_el_2_test: findpts_el_2_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
findpts_el_3_test2: rand_elt_test.o findpts_el_3_test2.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
findpts_el_3_test: findpts_el_3_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
findpts_local_test: rand_elt_test.o findpts_local_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
findpts_test: rand_elt_test.o findpts_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@

