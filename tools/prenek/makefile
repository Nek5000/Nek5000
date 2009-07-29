FLAGS =  $(FFLAGS_IN)
CFLAGS = $(CFLAGS_IN)
LFLAGS = $(LFLAGS_IN)
prefix = $(bin_nek_tools)

C     = code
L     = .
X     = /usr/X11R6/lib


NOBJS = prenek.o curve.o edit.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib.o subs.o zipper2.o postnek6.o screen.o revert.o crs.o mxm.o xdriver.o 

all: prex


prex:	$(NOBJS)
	$(F77) -o $(prefix)/prex $(NOBJS) -L$(X) -lX11 -lm

clean:
	'rm' *.o

plot.o		: plot.f 	basics.inc	; $(F77) -c $(FLAGS) plot.f
blas.o		: blas.f 	basics.inc	; $(F77) -c $(FLAGS) blas.f
screen.o	: screen.f	basics.inc	; $(F77) -c $(FLAGS) screen.f
mxm.o		: mxm.f		basics.inc	; $(F77) -c $(FLAGS) mxm.f
crs.o		: crs.f		basics.inc	; $(F77) -c $(FLAGS) crs.f
bound.o		: bound.f 	basics.inc	; $(F77) -c $(FLAGS) bound.f
prenek.o	: prenek.f	basics.inc	; $(F77) -c $(FLAGS) prenek.f
zipper.o	: zipper.f	basics.inc	; $(F77) -c $(FLAGS) zipper.f
zipper2.o	: zipper2.f	basics.inc	; $(F77) -c $(FLAGS) zipper2.f
curve.o		: curve.f 	basics.inc	; $(F77) -c $(FLAGS) curve.f
edit.o		: edit.f 	basics.inc	; $(F77) -c $(FLAGS) edit.f
build.o		: build.f 	basics.inc	; $(F77) -c $(FLAGS) build.f
xinterface.o: xinterface.f 	basics.inc	; $(F77) -c $(FLAGS) xinterface.f
postnek6.o	: postnek6.f 	basics.inc	; $(F77) -c $(FLAGS) postnek6.f
glomod.o	: glomod.f 	basics.inc	; $(F77) -c $(FLAGS) glomod.f
legend.o	: legend.f	basics.inc	; $(F77) -c $(FLAGS) legend.f
vprops.o	: vprops.f	basics.inc	; $(F77) -c $(FLAGS) vprops.f
iolib.o		: iolib.f	basics.inc	; $(F77) -c $(FLAGS) iolib.f
subs.o		: subs.f	basics.inc	; $(F77) -c $(FLAGS) subs.f
build1.o	: build1.f	basics.inc	; $(F77) -c $(FLAGS) build1.f
build2.o	: build2.f	basics.inc	; $(F77) -c $(FLAGS) build2.f
g3d.o		: g3d.f		basics.inc	; $(F77) -c $(FLAGS) g3d.f
xdriver.o	: xdriver.c			; $(CC) -c $(CFLAGS) xdriver.c
revert.o	: revert.c			; $(CC) -c $(CFLAGS) revert.c
