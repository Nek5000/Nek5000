PREFIX    = $(bin_nek_tools)

all: amg_hypre

amg_hypre: amg_hypre.o ./hypre/lib/libHYPRE.a 
	$(CC) -o $(PREFIX)/$@ $^ -lm $(LFLAGS)

%.o: %.c
	$(CC) -DHAVE_CONFIG_H -DHYPRE_TIMING $(CFLAGS) -I./hypre/include -c $<

%.a:
	@cd hypre; ./install

clean:
	@rm -rf amg_hypre *.o hypre/lib hypre/include hypre/hypre* hypre/*.tar*

.PHONY: clean
