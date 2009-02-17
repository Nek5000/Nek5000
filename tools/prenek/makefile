FLAGS =  $(FFLAGS_IN)
CFLAGS = $(CFLAGS_IN)
LFLAGS = $(LFLAGS_IN)
prefix = $(bin_nek_tools)

C     = code
L     = .
X     = /usr/X11R6/lib


GOBJS = prenek.o curve.o edit.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib_no_graph.o subs.o zipper2.o postnek6.o screen.o revert.o mxm.o xdriver.o error.o ivec.o queue.o stack.o bss_malloc.o adj_list.o pythag.o tqli.o smbwr.o sparse_matrix.o rsb_driver.o rsb.o

NOBJS = prenek.o curve.o edit.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib.o subs.o zipper2.o postnek6.o screen.o revert.o mxm.o xdriver.o error.o ivec.o queue.o stack.o bss_malloc.o adj_list.o pythag.o tqli.o smbwr.o sparse_matrix.o rsb_driver.o rsb.o blas.o


all: prex


prex_no_graph:	$(GOBJS)
	$(F77) -o $(prefix)/prex_no_graph $(GOBJS) -L$(X) -lX11 -lm
prex:	$(NOBJS)
	$(F77) -o $(prefix)/prex $(NOBJS) -L$(X) -lX11 -lm

clean:
	'rm' *.o

plot.o		: plot.f 	basics.inc	; $(F77) -c $(FLAGS) plot.f
blas.o		: blas.f 	basics.inc	; $(F77) -c $(FLAGS) blas.f
screen.o	 	: screen.f	basics.inc	; $(F77) -c $(FLAGS) screen.f
mxm.o			: mxm.f		basics.inc	; $(F77) -c $(FLAGS) mxm.f
bound.o		: bound.f 	basics.inc	; $(F77) -c $(FLAGS) bound.f

prenek.o		: prenek.f	basics.inc	; $(F77) -c $(FLAGS) prenek.f
zipper.o		: zipper.f	basics.inc	; $(F77) -c $(FLAGS) zipper.f
zipper2.o	: zipper2.f	basics.inc	; $(F77) -c $(FLAGS) zipper2.f
curve.o		: curve.f 	basics.inc	; $(F77) -c $(FLAGS) curve.f
edit.o		: edit.f 	basics.inc	; $(F77) -c $(FLAGS) edit.f
build.o		: build.f 	basics.inc	; $(F77) -c $(FLAGS) build.f
xinterface.o: xinterface.f 	basics.inc	; $(F77) -c $(FLAGS) xinterface.f
postnek6.o	: postnek6.f 	basics.inc	; $(F77) -c $(FLAGS) postnek6.f
glomod.o		: glomod.f 	basics.inc	; $(F77) -c $(FLAGS) glomod.f
legend.o		: legend.f	basics.inc	; $(F77) -c $(FLAGS) legend.f
vprops.o		: vprops.f	basics.inc	; $(F77) -c $(FLAGS) vprops.f
iolib.o		: iolib.f	basics.inc	; $(F77) -c $(FLAGS) iolib.f
iolib_no_graph.o		: iolib_no_graph.f	basics.inc	; $(F77) -c $(FLAGS) iolib_no_graph.f
subs.o		: subs.f	basics.inc	; $(F77) -c $(FLAGS) subs.f
build1.o		: build1.f	basics.inc	; $(F77) -c $(FLAGS) build1.f
build2.o		: build2.f	basics.inc	; $(F77) -c $(FLAGS) build2.f
g3d.o			: g3d.f		basics.inc	; $(F77) -c $(FLAGS) g3d.f
rsb.o			: rsb.F		basics.inc	; $(F77) -c $(FLAGS) rsb.F


xdriver.o	: xdriver.c			; $(CC) -c $(CFLAGS) xdriver.c
revert.o		: revert.c			; $(CC) -c $(CFLAGS) revert.c

bss_malloc.o: $C/bss_malloc.c	; $(CC) -c  $(CFLAGS) $C/bss_malloc.c
error.o		: $C/error.c		; $(CC) -c  $(CFLAGS) $C/error.c
ivec.o		: $C/ivec.c			; $(CC) -c  $(CFLAGS) $C/ivec.c
queue.o		: $C/queue.c		; $(CC) -c  $(CFLAGS) $C/queue.c
stack.o		: $C/stack.c		; $(CC) -c  $(CFLAGS) $C/stack.c

pythag.o		: $L/pythag.c		; $(CC) -c -I$C $(CFLAGS) $L/pythag.c
tqli.o		: $L/tqli.c			; $(CC) -c -I$C $(CFLAGS) $L/tqli.c
smbwr.o		: $L/smbwr.c		; $(CC) -c -I$C $(CFLAGS) $L/smbwr.c
adj_list.o	: $L/adj_list.c	; $(CC) -c -I$C $(CFLAGS) $L/adj_list.c
rsb_driver.o: $L/rsb_driver.c	; $(CC) -c -I$C $(CFLAGS) $L/rsb_driver.c
sparse_matrix.o	: $L/sparse_matrix.c	; $(CC) -c -I$C $(CFLAGS) $L/sparse_matrix.c
