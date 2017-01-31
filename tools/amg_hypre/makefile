CFLAGS    = -DHAVE_CONFIG_H -DHYPRE_TIMING 
LFLAGS    = $(LFLAGS_IN) -lm
PREFIX    = $(bin_nek_tools)

all: amg_hypre

amg_hypre: amg_hypre.o 
	$(CC) -o $(PREFIX)/$@ $^ $(LFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	-rm amg_hypre *.o

.PHONY: clean
