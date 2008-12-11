S     = .
L     = 
LIBS  = 
X     = /usr/X11R6/lib
P     = 

FLAGS =  
OLAGS = 
prefix = $(bin_nek_tools)


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
	'rm' postx


postnek.o  : postnek.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek.f 
postnek2.o : postnek2.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek2.f
postnek3.o : postnek3.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek3.f
locglob.o  : locglob.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  locglob.f 
postnek5.o : postnek5.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek5.f
postnek5a.o : postnek5a.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek5a.f
postnek6.o : postnek6.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek6.f
postnek7.o : postnek7.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek7.f
postnek8.o : postnek8.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek8.f
postnek9.o : postnek9.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  postnek9.f
speclib.o  : speclib.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  speclib.f
getfld.o   : getfld.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  getfld.f 
plot.o	   : plot.f 	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  plot.f
coef.o	   : coef.f 	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  coef.f
tekplot.o  : tekplot.f 	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  tekplot.f
legend.o   : legend.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  legend.f
animate.o  : animate.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  animate.f
subs.o     : subs.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  subs.f
g3d.o      : g3d.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  g3d.f
mxm.o      : mxm.f	basics.inc basicsp.inc	; $(F77) -c $P $(OLAGS)  mxm.f
genxyz.o   : genxyz.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  genxyz.f
screen.o   : screen.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  screen.f
scrdmp.o   : scrdmp.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  scrdmp.f
iolib.o    : iolib.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  iolib.f
iolib_no_graph.o    : iolib_no_graph.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  iolib_no_graph.f
tsort.o    : tsort.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  tsort.f
ssyev.o    : ssyev.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  ssyev.f
userf.o    : userf.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  userf.f
trap.o     : trap.f	basics.inc basicsp.inc	; $(F77) -c $P $(FLAGS)  trap.f
xinterface.o	: xinterface.f 	basics.inc	; $(F77) -c $P $(FLAGS)  xinterface.f
reverc.o	: revert.f			; $(F77) -c $P $(FLAGS)  revert.f
revert.o	: revert.c			; $(CC)  -c $(FLAGS)  revert.c
byte.o	: byte.c			; $(CC)  -c -O2   byte.c
#  xdriver.o	: xdriver.c		; $(CC)  -c -O2  -DLAPTOP xdriver.c
xdriver.o	: xdriver.c		; $(CC)  -c -O2   xdriver.c
#  xdriver.o	: xdriver.c		; cp  xdriver.ok xdriver.o
blas.o	   : ../../nek/blas.f	; $(F77) -c $P $(FLAGS)  ../../nek/blas.f
