S     = .
LIB   = .
P     = 
FLAGSF = -r8 -Mextend -Mlfs 
FLAGSC = 
FLAGS = -mcmodel=medium -O2
prefix = $(bin_nek_tools)

OBJS = main.o subs.o byte.o interp.o bound.o mxm.o speclib.o getfld.o dcopy.o

all: interp


interp:	$(OBJS)
	$(F77) $(FLAGS) -o $(prefix)/interp $(OBJS)  

clean:
	'rm' *.o

main.o        : main.f SIZE		 ;  $(F77) -c $P $(FLAGSF) main.f 
interp.o      : interp.f SIZE     ;  $(F77) -c $P $(FLAGSF) interp.f 
getfld.o      : getfld.f SIZE     ;  $(F77) -c $P $(FLAGSF) getfld.f 
bound.o       : bound.f  SIZE     ;  $(F77) -c $P $(FLAGSF) bound.f 
speclib.o     : speclib.f         ;  $(F77) -c $P $(FLAGSF) speclib.f 
subs.o        : subs.f            ;  $(F77) -c $P $(FLAGSF) subs.f 
mxm.o  	     : mxm.f             ;  $(F77) -c $P $(FLAGSF) mxm.f 
dcopy.o       : dcopy.f           ;  $(F77) -c $P $(FLAGSF) dcopy.f 
byte.o        : byte.c            ;  $(CC) -c     $(FLAGSC) byte.c
