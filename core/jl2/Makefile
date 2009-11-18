DATA_TYPES = -DGLOBAL_LONG
#JAMES_FLAGS=-DPREFIX=jl_ -DFPREFIX=jlf_ -DMPI -DUNDERSCORE $(DATA_TYPES)
#JAMES_FLAGS=-DPREFIX=jl_ -DMPI $(DATA_TYPES)
JAMES_FLAGS=-DPREFIX=jl_ $(DATA_TYPES)
SHELL=/bin/sh
CC=gcc
LIBS=-lm
GCC_WARN=-W -Wall -Wno-unused-function -Wno-unused-parameter
CEXTRA = $(GCC_WARN) -O3 -g $(JAMES_FLAGS)
#CEXTRA=-Minform=warn -O0 -g $(JAMES_FLAGS)
CCCMD=$(CC) $(CFLAGS) $(CEXTRA) -c
LINKCMD=$(CC) $(LDFLAGS) $(LIBS)
#RLINKCMD = $(LD) -r
.PHONY: cmds deps tests clean objects odepinfo

TESTS=sort_test sort_test2 sarray_sort_test spchol_test \
      comm_test crystal_test sarray_transfer_test \
      gs_test gs_test_old gs_unique_test \
      xxt_test xxt_test2
      
CRS=$(XXT)

tests: $(TESTS) ;
clean: ; @$(RM) $(TESTS) *.o *.s

cmds: ; @echo CC = $(CCCMD); echo LINK = $(LINKCMD);

deps: ; ./cdep.py *.c > makefile.cdep;

odepinfo: deps objects; @./odep_info.py *.o

-include makefile.cdep

%.o: %.c ; @echo CC $<; $(CCCMD) $< -o $@
objects: $(OBJECTS) ;

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
