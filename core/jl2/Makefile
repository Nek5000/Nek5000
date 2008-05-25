DATA_TYPES = -DGLOBAL_LONG
JAMES_FLAGS=-DPREFIX=jl_ -DMPI -DUNDERSCORE $(DATA_TYPES)
#JAMES_FLAGS=-DPREFIX=jl_ $(DATA_TYPES)
SHELL=/bin/sh
CC=mpicc
#CC = gcc-4.1
LIBS=-lm
GCC_WARN=-W -Wall -Wno-unused-function -Wno-unused-parameter
CEXTRA = $(GCC_WARN) -O0 -g $(JAMES_FLAGS)
#CEXTRA=-Minform=warn -O0 -g $(JAMES_FLAGS)
CCCMD=$(CC) $(CFLAGS) $(CEXTRA) -c
LINKCMD=$(CC) $(LDFLAGS) $(LIBS)
#RLINKCMD = $(LD) -r
.PHONY: default tests clean

TESTS=sort_test sarray_sort_test comm_test crystal_test sarray_transfer_test \
      gs_test gs_test2 gs_test_old gs_unique_test \
      spchol_test xxt_test xxt_test2

CRS=$(XXT)

default: tests

tests: $(TESTS)

clean:
	@$(RM) $(TESTS) *.o *.s

.c.o:
	$(CCCMD) $< -o $@

errmem.o: errmem.c
sort.o: sort.c sort_imp.c name.h errmem.h types.h
sort_test.o: sort_test.c name.h errmem.h types.h sort.h
comm.o: comm.c name.h errmem.h types.h comm.h
comm_test.o: comm_test.c name.h errmem.h types.h comm.h
crystal.o: crystal.c name.h errmem.h types.h comm.h
crystal_test.o: crystal_test.c name.h errmem.h types.h comm.h crystal.h
sarray_sort.o: sarray_sort.c name.h errmem.h types.h sort.h
sarray_sort_test.o: sarray_sort_test.c name.h errmem.h types.h sarray_sort.h
sarray_transfer.o: sarray_transfer.c name.h errmem.h types.h comm.h crystal.h sort.h
sarray_transfer_test.o: sarray_transfer_test.c name.h errmem.h types.h comm.h crystal.h sarray_transfer.h

gs_local.o: gs_local.c gs_defs.h name.h types.h
gs.o: gs.c gs_defs.h gs_local.h name.h errmem.h types.h comm.h crystal.h sarray_sort.h sarray_transfer.h
gs_test.o: gs_test.c name.h errmem.h types.h comm.h crystal.h gs.h gs_defs.h
gs_test2.o: gs_test2.c name.h errmem.h types.h comm.h crystal.h gs.h gs_defs.h
gs_test_old.o: gs_test_old.c types.h name.h
gs_unique_test.o: gs_unique_test.c name.h errmem.h types.h comm.h crystal.h gs.h gs_defs.h
GS_OBJECTS=gs.o gs_local.o \
sarray_sort.o sarray_transfer.o crystal.o comm.o sort.o errmem.o

sparse_cholesky.o: sparse_cholesky.c name.h sort.h errmem.h types.h
spchol_test.o: spchol_test.c name.h types.h errmem.h sparse_cholesky.h
xxt.o: xxt.c name.h types.h errmem.h sarray_sort.h gs.h comm.h gs_defs.h sparse_cholesky.h
XXT=sparse_cholesky.o xxt.o

fcrs.o: fcrs.c name.h types.h errmem.h comm.h crs.h
xxt_test.o: xxt_test.c name.h types.h errmem.h comm.h xxt.h
xxt_test2.o: xxt_test2.c name.h types.h errmem.h comm.h xxt.h

sort_test: sort_test.o sort.o errmem.o;					$(LINKCMD) $^ -o $@
sarray_sort_test: sarray_sort_test.o sarray_sort.o sort.o errmem.o;	$(LINKCMD) $^ -o $@
comm_test: comm_test.o comm.o errmem.o gs_local.o;			$(LINKCMD) $^ -o $@
crystal_test: crystal_test.o crystal.o errmem.o;			$(LINKCMD) $^ -o $@
sarray_transfer_test: sarray_transfer_test.o sarray_transfer.o crystal.o sort.o errmem.o;	$(LINKCMD) $^ -o $@
gs_test: gs_test.o $(GS_OBJECTS);					$(LINKCMD) $^ -o $@
gs_test2: gs_test2.o $(GS_OBJECTS);					$(LINKCMD) $^ -o $@
gs_test_old: gs_test_old.o $(GS_OBJECTS);				$(LINKCMD) $^ -o $@
gs_unique_test: gs_unique_test.o $(GS_OBJECTS);				$(LINKCMD) $^ -o $@
spchol_test: spchol_test.o sparse_cholesky.o errmem.o sort.o;		$(LINKCMD) $^ -o $@
xxt_test: xxt_test.o $(CRS) $(GS_OBJECTS);			$(LINKCMD) $^ -o $@
xxt_test2: xxt_test2.o $(CRS) $(GS_OBJECTS);			$(LINKCMD) $^ -o $@
