########################################################################
# Compiler and external dependences
########################################################################
CC        = mpicc
HYPRE_DIR = 'PathToHypreSource'/src/hypre

########################################################################
# Compiling and linking options
########################################################################
COPTS     = -O3 -Wall # -g
CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(CFLAGS_IN) $(COPTS) $(CINCLUDES) $(CDEFS)

LINKOPTS  = $(COPTS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++

########################################################################
# Rules for compiling the source files
########################################################################
all: amg_hypre

amg_hypre: amg_hypre.o 
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	-rm amg_hypre *.o

.PHONY: clean
