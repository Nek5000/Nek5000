S     = .
L     = .
LIBS  = ../../blas
C     = code
X     = /usr/X11R6/lib
F77   = pgf77 -Ktrap=fp
P     = 
CC    = pgcc

FLAGS = -g -O0 -Dr8

GOBJS = prenek.o edit.o curve.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib_no_graph.o subs.o zipper2.o postnek6.o screen.o revert.o mxm.o xdriver.o error.o ivec.o queue.o stack.o bss_malloc.o adj_list.o pythag.o tqli.o smbwr.o sparse_matrix.o rsb_driver.o rsb.o

NOBJS = prenek.o edit.o curve.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib.o subs.o zipper2.o postnek6.o screen.o revert.o mxm.o xdriver.o error.o ivec.o queue.o stack.o bss_malloc.o adj_list.o pythag.o tqli.o smbwr.o sparse_matrix.o rsb_driver.o rsb.o


all: prex print


prex_no_graph:	$(GOBJS)
	$(F77) -o prex_no_graph $(GOBJS) -L$(X) -lX11 -lm -L$(LIBS) -lblas
prex:	$(NOBJS)
	$(F77) -o prex $(NOBJS) -L$(X) -lX11 -lm -L$(LIBS) -lblas # -Bstatic

print:
	size prex


clean:
	'rm' *.o
	'rm' prex




plot.o		: plot.f 	basics.inc	; $(F77) $P -c $(FLAGS) plot.f
blas.o		: blas.f 	basics.inc	; $(F77) $P -c $(FLAGS) blas.f
screen.o	: screen.f	basics.inc	; $(F77) $P -c $(FLAGS) screen.f
mxm.o		: mxm.f		basics.inc	; $(F77) $P -c $(FLAGS) mxm.f
bound.o		: bound.f 	basics.inc	; $(F77) $P -c $(FLAGS) bound.f

prenek.o	: prenek.f	basics.inc	; $(F77) $P -c $(FLAGS) prenek.f
zipper.o	: zipper.f	basics.inc	; $(F77) $P -c $(FLAGS) zipper.f
zipper2.o	: zipper2.f	basics.inc	; $(F77) $P -c $(FLAGS) zipper2.f
edit.o		: edit.f 	basics.inc	; $(F77) $P -c $(FLAGS) edit.f
curve.o		: curve.f 	basics.inc	; $(F77) $P -c $(FLAGS) curve.f
build.o		: build.f 	basics.inc	; $(F77) $P -c $(FLAGS) build.f
xinterface.o	: xinterface.f 	basics.inc	; $(F77) $P -c $(FLAGS) xinterface.f
postnek6.o	: postnek6.f 	basics.inc	; $(F77) $P -c $(FLAGS) postnek6.f
glomod.o	: glomod.f 	basics.inc	; $(F77) $P -c $(FLAGS) glomod.f
legend.o	: legend.f	basics.inc	; $(F77) $P -c $(FLAGS) legend.f
vprops.o	: vprops.f	basics.inc	; $(F77) $P -c $(FLAGS) vprops.f
iolib.o		: iolib.f	basics.inc	; $(F77) $P -c $(FLAGS) iolib.f
iolib_no_graph.o		: iolib_no_graph.f	basics.inc	; $(F77) $P -c $(FLAGS) iolib_no_graph.f
subs.o		: subs.f	basics.inc	; $(F77) $P -c $(FLAGS) subs.f
build1.o	: build1.f	basics.inc	; $(F77) $P -c $(FLAGS) build1.f
build2.o	: build2.f	basics.inc	; $(F77) $P -c $(FLAGS) build2.f
g3d.o		: g3d.f		basics.inc	; $(F77) $P -c $(FLAGS) g3d.f
rsb.o		: rsb.F		basics.inc	; $(F77) $P -c $(FLAGS) rsb.F


xdriver.o	: xdriver.c		; $(CC) -c $(FLAGS) xdriver.c
revert.o	: revert.c		; $(CC) -c $(FLAGS) revert.c

bss_malloc.o	: $C/bss_malloc.c	; $(CC) -c  $(FLAGS) $C/bss_malloc.c
error.o		: $C/error.c		; $(CC) -c  $(FLAGS) $C/error.c
ivec.o		: $C/ivec.c		; $(CC) -c  $(FLAGS) $C/ivec.c
queue.o		: $C/queue.c		; $(CC) -c  $(FLAGS) $C/queue.c
stack.o		: $C/stack.c		; $(CC) -c  $(FLAGS) $C/stack.c

pythag.o	: $L/pythag.c		; $(CC) -c -I$C $(FLAGS) $L/pythag.c
tqli.o		: $L/tqli.c		; $(CC) -c -I$C $(FLAGS) $L/tqli.c
smbwr.o		: $L/smbwr.c		; $(CC) -c -I$C $(FLAGS) $L/smbwr.c
adj_list.o	: $L/adj_list.c		; $(CC) -c -I$C $(FLAGS) $L/adj_list.c
rsb_driver.o	: $L/rsb_driver.c	; $(CC) -c -I$C $(FLAGS) $L/rsb_driver.c
sparse_matrix.o	: $L/sparse_matrix.c	; $(CC) -c -I$C $(FLAGS) $L/sparse_matrix.c
