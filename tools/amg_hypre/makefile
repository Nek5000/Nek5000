PREFIX    = $(bin_nek_tools)

all: amg_hypre

amg_hypre: amg_hypre.o 
	$(CC) -o $(PREFIX)/$@ $^ -lm $(LFLAGS)

%.o: %.c
	$(CC) -DHAVE_CONFIG_H -DHYPRE_TIMING $(CFLAGS) -c $<

clean:
	-rm amg_hypre *.o

.PHONY: clean
