DATA_TYPES = -std=c99 -DGLOBAL_LONG_LONG
#JAMES_FLAGS=-DPREFIX=jl_ -DFPREFIX=jlf_ -DMPI -DUNDERSCORE $(DATA_TYPES)
JAMES_FLAGS=-DPREFIX=jl_ -DMPI $(DATA_TYPES)
#JAMES_FLAGS=-DPREFIX=jl_ $(DATA_TYPES)
SHELL=/bin/sh
CC=mpicc
LIBS=-lm
GCC_WARN=-W -Wall -Wno-unused-function -Wno-unused-parameter
CEXTRA = $(GCC_WARN) -O3 -g $(JAMES_FLAGS)
#CEXTRA=-Minform=warn -O0 -g $(JAMES_FLAGS)
#CBLAS=
#LBLAS=
CBLAS=-DUSE_CBLAS
LBLAS=-lcblas
CCCMD=$(CC) $(CFLAGS) $(CBLAS) $(CEXTRA)
LINKCMD=$(CC) $(LDFLAGS) $(LBLAS) $(LIBS)
#RLINKCMD = $(LD) -r
.PHONY: cmds deps tests clean objects odepinfo

TESTS=sort_test sort_test2 sarray_sort_test spchol_test \
      comm_test crystal_test sarray_transfer_test \
      gs_test gs_test_old gs_unique_test \
      xxt_test xxt_test2 \
      poly_test poly_test2 lob_bnd_test obbox_test

CRS=$(XXT)

tests: $(TESTS) ;
clean: ; @$(RM) $(TESTS) *.o *.s

cmds: ; @echo CC = $(CCCMD); echo LINK = $(LINKCMD);

deps: ; ./cdep.py *.c > makefile.cdep;

odepinfo: deps objects; @./odep_info.py *.o

-include makefile.cdep

%.o: %.c ; @echo CC $<; $(CCCMD) -c $< -o $@
%.s: %.c ; @echo CC -S $<; $(CCCMD) -S $< -o $@
objects: $(OBJECTS) ;

poly_imp.h: gen_poly_imp.c
	$(RM) poly_imp.h;
	$(CC) -lgmp -lm gen_poly_imp.c -o gen_poly_imp;
	./gen_poly_imp > poly_imp.h;
	$(RM) gen_poly_imp

GS_OBJECTS=gs.o sort.o sarray_transfer.o sarray_sort.o \
           gs_local.o fail.o crystal.o comm.o tensor.o

XXT=sparse_cholesky.o xxt.o

#AMG_DUMP=amg_dump.o
#AMG=amg.o

sort_test: sort.o fail.o sort_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
sort_test2: sort.o fail.o sort_test2.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
sarray_sort_test: sort.o fail.o sarray_sort.o sarray_sort_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
spchol_test: sparse_cholesky.o sort.o fail.o spchol_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
comm_test: comm.o tensor.o gs_local.o comm_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
crystal_test: fail.o crystal.o crystal_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
sarray_transfer_test: sarray_transfer.o sort.o fail.o crystal.o sarray_transfer_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@

gs_test: gs_test.o $(GS_OBJECTS);		@echo LINK $@; $(LINKCMD) $^ -o $@
gs_test_old: gs_test_old.o $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
gs_unique_test: gs_unique_test.o $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
xxt_test: xxt_test.o $(CRS) $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@
xxt_test2: xxt_test2.o $(CRS) $(GS_OBJECTS);	@echo LINK $@; $(LINKCMD) $^ -o $@

poly_test: poly.o fail.o poly_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
poly_test2: poly.o fail.o poly_test2.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
lob_bnd_test: tensor.o poly.o lob_bnd.o fail.o lob_bnd_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
obbox_test: poly.o obbox.o tensor.o lob_bnd.o fail.o obbox_test.o ; @echo LINK $@; $(LINKCMD) $^ -o $@
