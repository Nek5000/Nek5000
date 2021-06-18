### makefile user
OBJDIR=obj
S:=/Users/yiminllin/Desktop/github/Nek5000/
CC=mpicc
CFLAGS+=
PPS=NOMPIIO MPI UNDERSCORE GLOBAL_LONG_LONG TIMER
PPS_C = $(patsubst %,-D%,$(PPS))
L2 = $(G) -O2
OPT_INCDIR=./
cFL2   = $(L2) $(CFLAGS) $(PPS_C) -I$S/core -I$(OPT_INCDIR)

nekio.o:$S/bin/nekio.c;     $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $< -o $@

