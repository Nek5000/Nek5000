prefix = $(bin_nek_tools)

X = /usr/lib/X11
UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
 X = /usr/lib/X11
else 
ifeq ($(UNAME),Darwin)
  VERSION1 := $(shell sw_vers -productVersion)
  ifeq "10.8" "$(word 1, $(sort 10.8 $(VERSION1)))"
     CFLAGS += -I/opt/X11/include
     X      = /opt/X11/lib
  else
     CFLAGS += -I/usr/X11/include
     X      = /usr/X11/lib
  endif
endif
endif

L     = .


NOBJS = prenek.o curve.o edit.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib.o subs.o zipper2.o postnek6.o screen.o revert.o crs.o mxm.o xdriver.o 

TOBJS = prenek.o curve.o edit.o build.o build1.o build2.o bound.o plot.o xinterface.o glomod.o legend.o vprops.o iolib_no_graph.o subs.o zipper2.o postnek6.o screen.o revert.o crs.o mxm.o xdriver.o


all: prex pretex

prex:	$(NOBJS)
	$(FC) -o $(prefix)/prex $(NOBJS) -L$(X) -lX11 -lm
	@echo $(UNAME)

pretex:	$(TOBJS)
	$(FC) -o $(prefix)/pretex $(TOBJS) -L$(X) -lX11 -lm

clean:
	@rm -f *.o
	
plot.o		: plot.f 	basics.inc	; $(FC) -c $(FFLAGS) plot.f
iolib_no_graph.o: iolib_no_graph.f basics.inc   ; $(FC) -c $(FFLAGS) iolib_no_graph.f
blas.o		: blas.f 	basics.inc	; $(FC) -c $(FFLAGS) blas.f
screen.o	: screen.f	basics.inc	; $(FC) -c $(FFLAGS) screen.f
mxm.o		: mxm.f		basics.inc	; $(FC) -c $(FFLAGS) mxm.f
crs.o		: crs.f		basics.inc	; $(FC) -c $(FFLAGS) crs.f
bound.o		: bound.f 	basics.inc	; $(FC) -c $(FFLAGS) bound.f
prenek.o	: prenek.f	basics.inc	; $(FC) -c $(FFLAGS) prenek.f
zipper.o	: zipper.f	basics.inc	; $(FC) -c $(FFLAGS) zipper.f
zipper2.o	: zipper2.f	basics.inc	; $(FC) -c $(FFLAGS) zipper2.f
curve.o		: curve.f 	basics.inc	; $(FC) -c $(FFLAGS) curve.f
edit.o		: edit.f 	basics.inc	; $(FC) -c $(FFLAGS) edit.f
build.o		: build.f 	basics.inc	; $(FC) -c $(FFLAGS) build.f
xinterface.o    : xinterface.f 	basics.inc	; $(FC) -c $(FFLAGS) xinterface.f
postnek6.o	: postnek6.f 	basics.inc	; $(FC) -c $(FFLAGS) postnek6.f
glomod.o	: glomod.f 	basics.inc	; $(FC) -c $(FFLAGS) glomod.f
legend.o	: legend.f	basics.inc	; $(FC) -c $(FFLAGS) legend.f
vprops.o	: vprops.f	basics.inc	; $(FC) -c $(FFLAGS) vprops.f
iolib.o		: iolib.f	basics.inc	; $(FC) -c $(FFLAGS) iolib.f
subs.o		: subs.f	basics.inc	; $(FC) -c $(FFLAGS) subs.f
build1.o	: build1.f	basics.inc	; $(FC) -c $(FFLAGS) build1.f
build2.o	: build2.f	basics.inc	; $(FC) -c $(FFLAGS) build2.f
g3d.o		: g3d.f		basics.inc	; $(FC) -c $(FFLAGS) g3d.f
xdriver.o	: xdriver.c			; $(CC) -c $(CFLAGS) xdriver.c
revert.o	: revert.c			; $(CC) -c $(CFLAGS) revert.c
