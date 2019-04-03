PREFIX = $(bin_nek_tools)
LIBS = ./hypre/lib/libHYPRE.a

all: lib nekamg_setup

nekamg_setup: amg_setup.o
	$(CC) -o $(PREFIX)/$@ $^ $(LIBS) $(LDFLAGS) -lm

%.o: %.c
	$(CC) -DHAVE_CONFIG_H -DHYPRE_SEQUENTIAL -DHYPRE_TIMING $(CFLAGS) -I./hypre/include -c $<

lib:
	@cd hypre; env CC="$(CC)" CFLAGS="$(BIGMEM)" FC="$(FC)" FFLAGS="$(BIGMEM)" ./install

clean:
	@rm -rf amg_setup *.o hypre/lib hypre/include hypre/hypre*

.PHONY: clean
