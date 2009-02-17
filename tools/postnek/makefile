FLAGS =  $(FFLAGS_IN)
CFLAGS = $(CFLAGS_IN)
LFLAGS = $(LFLAGS_IN)
prefix = $(bin_nek_tools)

S     = .
X     = /usr/X11R6/lib

NOBJS = postnek.o postnek2.o postnek3.o postnek5.o postnek6.o tsort.o \
	postnek8.o postnek9.o plot.o getfld.o legend.o userf.o revert.o \
	trap.o animate.o genxyz.o screen.o g3d.o subs.o xinterface.o \
	locglob.o postnek5a.o blas.o \
	xdriver.o scrdmp.o coef.o postnek7.o speclib.o mxm.o byte.o ssyev.o

all: postx


## posty:	$(NOBJS) iolib_no_graph.o
## 	$(F77) -o  $(prefix)/posty $(NOBJS) iolib_no_graph.o -L$(X) -lX11 -lm 

postx:	$(NOBJS) iolib.o
	$(F77) -o  $(prefix)/postx $(NOBJS) iolib.o -L$(X) -lX11 -lm

print:
	size postx

clean:
	'rm' *.o

postnek.o  : postnek.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek.f 
postnek2.o : postnek2.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek2.f
postnek3.o : postnek3.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek3.f
locglob.o  : locglob.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  locglob.f 
postnek5.o : postnek5.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek5.f
postnek5a.o : postnek5a.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek5a.f
postnek6.o : postnek6.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek6.f
postnek7.o : postnek7.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek7.f
postnek8.o : postnek8.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek8.f
postnek9.o : postnek9.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  postnek9.f
speclib.o  : speclib.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  speclib.f
getfld.o   : getfld.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  getfld.f 
plot.o	   : plot.f 	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  plot.f
coef.o	   : coef.f 	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  coef.f
tekplot.o  : tekplot.f 	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  tekplot.f
legend.o   : legend.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  legend.f
animate.o  : animate.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  animate.f
subs.o     : subs.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  subs.f
g3d.o      : g3d.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  g3d.f
mxm.o      : mxm.f	basics.inc basicsp.inc	; $(F77) -c $(OLAGS)  mxm.f
genxyz.o   : genxyz.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  genxyz.f
screen.o   : screen.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  screen.f
scrdmp.o   : scrdmp.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  scrdmp.f
iolib.o    : iolib.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  iolib.f
iolib_no_graph.o    : iolib_no_graph.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  iolib_no_graph.f
tsort.o    : tsort.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  tsort.f
ssyev.o    : ssyev.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  ssyev.f
userf.o    : userf.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  userf.f
trap.o     : trap.f	basics.inc basicsp.inc	; $(F77) -c $(FLAGS)  trap.f
xinterface.o	: xinterface.f 	basics.inc	; $(F77) -c $(FLAGS)  xinterface.f
reverc.o	: revert.f			; $(F77) -c $(FLAGS)  revert.f
blas.o	   : ../../nek/blas.f	; $(F77) -c $(FLAGS)  ../../nek/blas.f

revert.o	: revert.c				; $(CC)  -c $(CFLAGS)  revert.c
byte.o   : ../../nek/byte.c   ; $(CC)  -c $(CFLAGS)  ../../nek/byte.c
#  xdriver.o	: xdriver.c		; $(CC)  -c -O2  -DLAPTOP xdriver.c
xdriver.o	: xdriver.c			; $(CC)  -c -O2   xdriver.c
#  xdriver.o	: xdriver.c		; cp  xdriver.ok xdriver.o
