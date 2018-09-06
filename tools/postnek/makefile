prefix = $(bin_nek_tools)

X = /usr/lib/X11
UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
 VERSION1 := $(shell uname -r)
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

S     = .

NOBJS = postnek.o postnek2.o postnek3.o postnek5.o postnek6.o tsort.o \
	postnek8.o postnek9.o plot.o getfld.o legend.o userf.o revert.o \
	trap.o animate.o genxyz.o screen.o g3d.o subs.o xinterface.o \
	locglob.o postnek5a.o blas.o \
	xdriver.o scrdmp.o coef.o postnek7.o speclib.o mxm.o byte.o ssyev.o \
	finiparser.o iniparser.o dictionary.o

all: postx postex

postx:	$(NOBJS) iolib.o
	$(FC) -o  $(prefix)/postx $(NOBJS) iolib.o -L$(X) -lX11 -lm
	@echo $(UNAME)
	@echo $(VERSION1)

postex:	$(NOBJS) iolib_no_graph.o
	$(FC) -o  $(prefix)/postex $(NOBJS) iolib_no_graph.o -L$(X) -lX11 -lm


print:
	size postx

clean:
	@rm -f *.o

postnek.o  : postnek.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek.f 
postnek2.o : postnek2.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek2.f
postnek3.o : postnek3.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek3.f
locglob.o  : locglob.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  locglob.f 
postnek5.o : postnek5.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek5.f
postnek5a.o : postnek5a.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek5a.f
postnek6.o : postnek6.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek6.f
postnek7.o : postnek7.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek7.f
postnek8.o : postnek8.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek8.f
postnek9.o : postnek9.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  postnek9.f
speclib.o  : speclib.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  speclib.f
getfld.o   : getfld.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  getfld.f 
plot.o	   : plot.f 	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  plot.f
coef.o	   : coef.f 	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  coef.f
tekplot.o  : tekplot.f 	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  tekplot.f
legend.o   : legend.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  legend.f
animate.o  : animate.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  animate.f
subs.o     : subs.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  subs.f
g3d.o      : g3d.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  g3d.f
mxm.o      : mxm.f	basics.inc basicsp.inc	; $(FC) -c $(OLAGS)  mxm.f
genxyz.o   : genxyz.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  genxyz.f
screen.o   : screen.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  screen.f
scrdmp.o   : scrdmp.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  scrdmp.f
iolib.o    : iolib.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  iolib.f
iolib_no_graph.o    : iolib_no_graph.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  iolib_no_graph.f
tsort.o    : tsort.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  tsort.f
ssyev.o    : ssyev.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  ssyev.f
userf.o    : userf.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  userf.f
trap.o     : trap.f	basics.inc basicsp.inc	; $(FC) -c $(FFLAGS)  trap.f
xinterface.o	: xinterface.f 	basics.inc	; $(FC) -c $(FFLAGS)  xinterface.f
reverc.o	: revert.f			; $(FC) -c $(FFLAGS)  revert.f
blas.o	   : ../../core/3rd_party/blas.f	; $(FC) -c $(FFLAGS)  ../../core/3rd_party/blas.f

finiparser.o	   : ../../core/3rd_party/finiparser.c	; $(CC) -c $(CFLAGS)  ../../core/3rd_party/finiparser.c
iniparser.o	   : ../../core/3rd_party/iniparser.c	; $(CC) -c $(CFLAGS)  ../../core/3rd_party/iniparser.c
dictionary.o	   : ../../core/3rd_party/dictionary.c	; $(CC) -c $(CFLAGS)  ../../core/3rd_party/dictionary.c

revert.o	: revert.c			; $(CC)  -c $(CFLAGS)  revert.c
byte.o   	: byte.c   			; $(CC)  -c $(CFLAGS)  byte.c
xdriver.o	: xdriver.c			; $(CC)  -c -O2  $(CFLAGS) xdriver.c
