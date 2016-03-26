# This include file is used by the makenek script
# to automatically create a makefile for Nek5000 
# (c) 2008,2009,2010 UCHICAGO ARGONNE, LLC
# ------------------------------------------------

echo "makenek - automatic build tool for Nek5000"

if [ "$PPLIST" == "?" ]; then
  echo "available pre-processor symbols:" 
  echo "  NOTIMER   disable runtime statistics"
  echo "  MPITIMER  enable MPI runtime statistics"
  echo "  MPIIO     use MPI-IO I/O kernel (experimental)"
  echo "  BG        enable Blue Gene optimizations (BG/L and BG/P)"
  echo "  K10_MXM   use optimized MxM kernel for AMD Family 10h processors" 
  echo "  CVODE     use ODE solver from Sundials to solve for IFIELD>1 (experimental)"
  echo "  MOAB      enable MOAB/CUBIT support (experimental)"
  echo "  NEKNEK    enable multidomain support (experimental)"
  echo "  NO_LAPACK disable nek-supplied LAPACK routines"
  echo "  BLAS_MXM  disable nek-supplied BLAS routines"
  echo "  EXTBAR    adds underscore to exit call(for BGQ)"
  exit 1
fi

# first do some checks ...
if [ $# -eq 0 ]; then
  echo ""
  echo "usage: makenek [.usr filename | clean] < -nocompile >"
  echo ""
  exit 1
fi

if [ "$mver" == "" ]; then
  echo "ERROR: Your makenek is too old!"
  echo "Please upgrade to last version shipped with the Nek source."
  exit 1
else
  if [ $mver -lt 1 ]; then
    echo "ERROR: Your makenek is too old!"
    echo "Please upgrade to last version shipped with the Nek source."
    exit 1
  fi
fi

if [ $1 = "clean" ]; then
  make clean
  exit 0
fi

NOCOMPILE=0
if [ "$2" == "-nocompile" ]; then
  NOCOMPILE=1
fi 

CASENAME=$1
CASEDIR=`pwd`
APATH_SRC=`cd $SOURCE_ROOT; pwd`
SOURCE_ROOT=$APATH_SRC
APATH_SRC=`cd $SOURCE_ROOT/cmt; pwd`
SOURCE_ROOT_CMT=$APATH_SRC
echo "CMT root" $SOURCE_ROOT_CMT

IFMOAB=false
IFNEKNEK=false

# do some basic checks
if [ "$CASEDIR" == "$SOURCE_ROOT" ]; then
   echo "FATAL ERROR: Working directory has to be different from the source!"
   exit 1
fi

if [ ! -f ${CASENAME}'.usr' ]; then
   echo "FATAL ERROR: Cannot find $1.usr!"
   exit 1
fi

if [ -f SIZEu ]; then
   if [ -f SIZE ]; then
      echo "FATAL ERROR: found SIZEu and SIZE!"
      exit 1
   fi
   NEK_WARN="WARNING: SIZEu will be deprecated in a future release. Please rename to SIZE!"
   export NEK_WARN 
   mv -v SIZEu SIZE
fi

if [ ! -f SIZE ]; then
   echo "FATAL ERROR: Cannot find SIZE!"
   exit 1
fi

if [ ! -f $SOURCE_ROOT/makefile.template ]; then
  echo "FATAL ERROR: Cannot find $SOURCE_ROOT/makefile.template!"
  exit 1
fi

# test F77 compiler
which `echo $F77 | awk '{print $1}'` 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $F77!"
  exit 1
fi
\rm test_f77.o 2>/dev/null
# basic compiler test
cat > test_f77.f << _ACEOF
      subroutine test
      end
_ACEOF
$F77 -c test_f77.f >/dev/null
if [ ! -f test_f77.o ]; then
  echo "FATAL ERROR: Basic compiler test for $F77 failed!"
  exit 1
fi
\rm test_f77.* 2>/dev/null

# test C compiler
which `echo $CC | awk '{print $1}'` 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $CC!"
  exit 1
fi
\rm test_cc.o 2>/dev/null
cat > test_cc.c << _ACEOF
      void function(){}
_ACEOF
$CC -c test_cc.c 1>/dev/null
if [ ! -f test_cc.o ]; then
  echo "FATAL ERROR: Basic compiler test for $CC failed!"
  exit 1
fi
\rm test_cc.* 2>/dev/null

# initial clean-up
rm -f nek5000 2>/dev/null
rm -f ./obj/subuser.o        2>/dev/null

# Check ptr size
cat > tmp.c << _ACEOF
#include <stdlib.h>
#include <stdio.h>
int main()
{
  int *p;printf("%li\n",sizeof(p));
}
_ACEOF
$CC tmp.c 2>&1>/dev/null
ptrSize=`./a.out`
rm tmp.c a.out
if [ "$ptrSize" == "8" ]
then
  PPLIST="${PPLIST} PTRSIZE8"
fi

# Check if the compiler adds an underscore to external functions
UNDERSCORE=false
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call byte_rewind
      end
_ACEOF
$F77 -c test_underscore.f 2>&1 >/dev/null 
nm test_underscore.o | grep byte_rewind_ 1>/dev/null
if [ $? -eq 0 ] 
then
  UNDERSCORE=true
fi
\rm test_underscore.* 2>/dev/null

F77comp=$F77
# mpi-wrappers?
echo $F77comp | grep -i mp >/dev/null
if [ $? -eq 0 ]; then   
# trying to figure out which compiler we're using 
  F77ok=0
  F77comp_=`$F77 -showme 2>/dev/null 1>.tmp` 
  F77comp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}'`
  if [ -f "`which $F77comp 2>/dev/null`" ]; then
    F77ok=1
  else
    F77comp_=`$F77 -show 2>/dev/null 1>.tmp`
    F77comp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}'`
    if [ -f "`which $F77comp 2>/dev/null`" ]; then
      F77ok=1
    fi
  fi
else
  F77comp=`echo $F77 | awk '{print $1}'`
  if [ -f "`which $F77comp 2>/dev/null`" ]; then
    F77ok=1
  fi
fi
\rm -f .tmp
if [ $F77ok -eq 0 ]; then
  F77comp="unknown"
fi

# assign F77 compiler specific flags
case $F77comp in
  *pgf*)        P="-r8 -Mpreprocess"
               ;;
  *gfortran*)   P="-fdefault-real-8 -fdefault-double-8 -x f77-cpp-input"
               ;;
  *ifort*)      P="-r8 -fpconstant -fpp"
               ;;
  *pathf*)      P="-r8 -cpp -fno-second-underscore"
               ;;
  *xlf*)       P="-qrealsize=8 -qdpc=e -qsuffix=cpp=f"
               PPPO="-WF,"
               F77="${F77} -qsuppress=cmpmsg"
               ;;
  *ftn*)        P="-r8 -Mpreprocess"
               ;;
  *sunf*)       P="-r8const -xtypemap=real:64 -fpp"
               ;;
  *open*)       P="-r8 -cpp -fno-second-underscore"
               ;;
  *)  echo "ERROR: Unable to detect compiler!"
      echo "        - don't know how to promote datatype REAL to 8 bytes"
      echo "        - don't know how to invoke the C pre-processor (CPP) before compilation"
      echo "       Please edit the makefile and specify the requested compiler flags using the P variable."
      echo ""
      P="<specify your compiler flags here>"
      NOCOMPILE=1
      read;;
esac
export PPPO

# set preprocessor symbols
if [ "$IFMPI" == "false" -o "$IFMPI" == "no" ]; then
  IFMPI=false
else
  # default
  IFMPI=true
  PPLIST="${PPLIST} MPI"
fi
export IFMPI

# Check size of long int
cat > tmp.c << _ACEOF
#include <stdlib.h>
#include <stdio.h>
int main()
{
  int i;
  i=sizeof(long int);
  printf("%i\n",i);
}
_ACEOF
$CC tmp.c 2>&1>/dev/null
longIntTest=`./a.out`
rm tmp.c a.out
if [ "$longIntTest" == "8" ]
then
  PPLIST="${PPLIST} LONGINT8"
fi

if [ "$UNDERSCORE" == "true" ]; then
  PPLIST="${PPLIST} UNDERSCORE"
fi 

if [ "$IFTIMING" == "false" -o "$IFTIMING" == "no" ]; then
    PPLIST="${PPLIST} NOTIMER"
fi

PPLIST="${PPLIST} GLOBAL_LONG_LONG"

# add user defined preprocessor symbols
for i in $PPLIST
do
   if [ "$i" == "MOAB" ]; then
     IFMOAB=true
     echo $F77comp | grep 'gfort' >/dev/null 
     if [ $? -eq 0 ]; then
       P="${P} -fcray-pointer"
     fi
     if [ "$MOAB_DIR" == "" ]; then
       echo "ABORT: Please specify path to MOAB in MOAB_DIR!" 
       exit 1
     else
       if [ ! -d $MOAB_DIR ]; then
         echo "ABORT: Cannot find " $MOAB_DIR 
         exit 1
       fi 
       if [ ! -f $MOAB_DIR/lib/iMesh-Defs.inc ]; then
         echo "ABORT: Cannot find iMesh-Defs.inc in" $MOAB_DIR"/lib; make sure MOAB is installed."
	 echo "       (to install, use 'make install' in" $MOAB_DIR")"
         exit 1
       fi 
     fi
   fi

if [ "$i" == "NEKNEK" ]; then
     IFNEKNEK=true
fi

   if [ "$i" == "CVODE" ]; then
    if [ "$CVODE_DIR" == "" ]; then
       echo "ABORT: Please specify path to CVODE in CVODE_DIR!" 
       exit 1
    else
       if [ ! -d $CVODE_DIR ]; then
         echo "ABORT: Cannot find " $CVODE_DIR 
         exit 1
       fi
    fi
    USR_LFLAGS="${USR_LFLAGS} -L$CVODE_DIR -lsundials_fcvode -lsundials_cvode"
    if [ "$IFMPI" == "true"  ]; then
      USR_LFLAGS="${USR_LFLAGS} -lsundials_fnvecparallel -lsundials_nvecparallel"
    else
        USR_LFLAGS="${USR_LFLAGS} -lsundials_fnvecserial -lsundials_nvecserial"
    fi
   fi
done

# AMG setup
if [ "$IFAMG_DUMP" == "true" ]; then
    if [ "$IFAMG" != "true" -a "$IFAMG" != "yes" ]; then
      echo "ABORT: Cannot dump AMG files if AMG is disabled!"
      exit 1
    fi
    IFAMG_DUMP=true
else
    IFAMG_DUMP=false
fi
export IFAMG_DUMP

if [ "$IFAMG" == "true" -o "$IFAMG" == "yes" ]; then
  IFAMG=true
else
  # default
  IFAMG=false
fi

# CMT flag
if [ "$IFCMT" == "true" -o "$IFCMT" == "yes" ]; then
  IFCMT=true
else
  # default
  IFCMT=false
fi
export IFCMT

# MXM-handling
if [ "$MXM_USER" != "" ]; then
    echo "ABORT: makenek version is too old! Update with latest from source!"
    exit 1 
fi
MXM_USER="mxm_std.o blas.o"
echo $PPLIST | grep 'BG' >/dev/null 
if [ $? -eq 0 ]; then
   MXM_USER="mxm_std.o bg_aligned3.o bg_mxm44.o bg_mxm44_uneven.o bg_mxm3.o blas.o" 
   OPT_FLAGS_STD="-qarch=450 -qtune=450"
   OPT_FLAGS_MAG="-O5 -qarch=450d -qtune=450"
fi
echo $PPLIST | grep 'K10_MXM' >/dev/null 
if [ $? -eq 0 ]; then
   MXM_USER="mxm_std.o k10_mxm.o blas.o" 
   USR_LFLAGS="${USR_LFLAGS} ${SOURCE_ROOT}/libk10_mxm.a"
fi
echo $PPLIST | grep 'BLAS_MXM' >/dev/null 
if [ $? -eq 0 ]; then
   MXM_USER="mxm_std.o" 
fi
LPACK="dsygv.o ssygv.o"
echo $PPLIST | grep 'NO_LAPACK' >/dev/null 
if [ $? -eq 0 ]; then
   LPACK=" "
fi

# set optimization flags
L0="\$(G) -O0"
L2="\$(G) -O2"
L3="\$(G) -O3" 
L4="\$(L3)"

# user specified opt flags
if [ "$OPT_FLAGS_STD" != "" ]; then
  echo $OPT_FLAGS_STD | grep "\-O." 1>/dev/null
  if [ $? -eq 0 ]; then
    L2="\$(G) $OPT_FLAGS_STD"
    L3="\$(G) $OPT_FLAGS_STD" 
  else
    L2="\$(G) -O2 $OPT_FLAGS_STD"
    L3="\$(G) -O3 $OPT_FLAGS_STD"
  fi
fi

if [ "$OPT_FLAGS_MAG" != "" ]; then
    L4="\$(G) $OPT_FLAGS_MAG"
fi

if [ "$USR_LIB" != "" ]; then
    USR_LFLAGS="${USR_LFLAGS} ${USR_LIB}"
fi

# tweak makefile template 
echo "generating makefile ..."
rm -rf makefile 2>/dev/null

sed -e "s:^F77[ ]*=.*:F77\:=$F77:" \
-e "s:^CC[ ]*=.*:CC\:=$CC:" \
-e "s:^G[ ]*=.*:G\:=$G:" \
-e "s:^OPT_FLAGS[ ]*=.*:OPT_FLAGS\:=$OPT_FLAGS:" \
-e "s/^P[ ]*=.*/P:=$P/" \
-e "s/^L0[ ]*=.*/L0=$L0/" \
-e "s/^L2[ ]*=.*/L2=$L2/" \
-e "s/^L3[ ]*=.*/L3=$L3/" \
-e "s/^L4[ ]*=.*/L4=$L4/" \
-e "s/^PPPO[ ]*=.*/PPPO=$PPPO/" \
-e "s/^PPS[ ]*=.*/PPS=$PPLIST/" \
-e "s:^MXM[ ]*=.*:MXM=$MXM_USER:" \
-e "s:^LPACK[ ]*=.*:LPACK=$LPACK:" \
-e "s/^IFAMG[ ]*=.*/IFAMG:=$IFAMG/" \
-e "s/^IFAMG_DUMP[ ]*=.*/IFAMG_DUMP:=$IFAMG_DUMP/" \
-e "s/^IFMPI[ ]*=.*/IFMPI:=$IFMPI/" \
-e "s/^IFCMT[ ]*=.*/IFCMT:=$IFCMT/" \
-e "s/^IFMOAB[ ]*=.*/IFMOAB:=$IFMOAB/" \
-e "s:^MOAB_DIR[ ]*=.*:MOAB_DIR\:=${MOAB_DIR}:" \
-e "s/^IFNEKNEK[ ]*=.*/IFNEKNEK:=$IFNEKNEK/" \
-e "s:^USR[ ]*=.*:USR\:=$USR:" \
-e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS\:=$USR_LFLAGS:" \
-e "s:^S[ ]*=.*:S\:=${SOURCE_ROOT}:" \
-e "s:^S2[ ]*=.*:S2\:=${SOURCE_ROOT_CMT}:" $SOURCE_ROOT/makefile.template >.makefile

echo $G | grep '\-g' 1>/dev/null
if [ $? -eq 0 ]; then
  sed 's/-O[1-4]/-O0/g' .makefile > .makefile.tmp
  mv .makefile.tmp .makefile
  echo "Activate DEBUG mode"
fi

if [ "$IFVISIT" == "true" -o "$IFVISIT" == "yes" ]; then
  echo ""
  echo "Activate VisIt Sim Mode"
  IFVISIT=true

  if [ -z "$VISIT_INC" ]; then
    VISIT_INC=${VISIT_INSTALL}/libsim/V2/include
  fi
  echo "  VISIT_INC = ($VISIT_INC)"
  if [ -z "$VISIT_LIB" ]; then
    VISIT_LIB=${VISIT_INSTALL}/libsim/V2/lib
  fi
  echo "  VISIT_LIB = ($VISIT_LIB)"
  if [ -n "$VISIT_STOP" ]; then
    VISIT_STOP="VISIT_STOP"
    echo "  VISIT_STOP is defined. Nek5000 will wait for VisIt to Connect."
  fi

  sed -e "s/^IFVISIT[ ]*=.*/IFVISIT:=${IFVISIT}/" \
      -e "s:^VISIT_IFLAGS[ ]*=.*:VISIT_IFLAGS\:=-I${VISIT_INC}:" \
      -e "s:^VISIT_LFLAGS[ ]*=.*:VISIT_LFLAGS\:=-L${VISIT_LIB} -lsimV2 -lsimV2f -ldl -lstdc++:" \
      -e "s/^PPS[ ]*=.*/PPS=${PPLIST} VISIT ${VISIT_STOP}/" \
      .makefile > .makefile.tmp

  mv .makefile.tmp .makefile
  echo ""
fi

if [ "$USR" != "" ]; then
  echo "###########################################################" >> makefile
  echo "include makefile_usr.inc" >> .makefile
fi

if [ -f .makefile ]; then
  sed -e "1i\\
### makefile automatically created by makenek `date +"%m/%d/%Y %T"` ###" \
-e "s:^CASEDIR[ ]*=.*:CASEDIR\:=${CASEDIR}:" \
-e "s:^CASENAME[ ]*=.*:CASENAME\:=${CASENAME}:" .makefile > makefile 
else
  echo "ERROR: Nek Makefile could not be created!"
  exit 1 
fi
\rm .makefile 2>/dev/null

# add plugins
for i in $PLUGIN_LIST
do
  if [ ! -f makenek_plugin_$i ]; then
     echo "ABORT: cannot find plug-in config file" makenek_plugin_$i
     exit 1
  fi
  source makenek_plugin_$i

  plugindir="SRC_PLUGIN"
  plugindir=`echo ${!plugindir}`

  if [ ! -f $plugindir/makenek.inc ]; then
     echo "ABORT: Cannot find plug-in include file file" $plugindir/makenek.inc
     exit 1
  fi
  source $plugindir/makenek.inc
  if [ $ierr -ne 0 ]; then
     echo "ABORT: plug-in config failed!"
     exit 1
  else
    echo "found plug-in" $i in $plugindir
  fi
done

# tweak SIZE file
if [ -f "./SIZE" ]; then
  cat SIZE | grep -i 'lxo' >/dev/null
  if [ $? -ne 0 ]; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)' >>SIZE
  fi
  cat SIZE | grep -i 'lpart' >/dev/null
  if [ $? -ne 0 ]; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter(lpart = 1  ) ! max number of particles/proc' >>SIZE
  fi
  cat SIZE | grep -i 'ax1' >/dev/null
  if [ $? -ne 0 ]; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      integer ax1,ay1,az1,ax2,ay2,az2' >> SIZE
     echo '      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! running averages' >> SIZE
  fi
  cat SIZE | grep -i 'lys=lxs' >/dev/null
  if [ $? -ne 0 ]; then
     cat SIZE | grep -iv lxs > SIZE.x; mv SIZE.x SIZE  # Clean existing SIZE file of old version
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Preconditioner' >> SIZE

  fi
  cat SIZE | grep -i 'lfdm' >/dev/null
  if [ $? -ne 0 ]; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (lfdm=0)  ! == 1 for fast diagonalization method' >> SIZE
  fi
  cat SIZE | grep -i 'nio' >/dev/null
  if [ $? -ne 0 ]; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      common/IOFLAG/nio  ! for logfile verbosity control' >> SIZE
  fi

  echo $PPLIST | grep 'CVODE' >/dev/null
  if [ $? -eq 0 ]; then
  cat SIZE | grep -i 'cv_maxl' >/dev/null
  if [ $? -ne 0 ]; then
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      integer cv_maxl' >> SIZE
      echo '      parameter(cv_maxl = 20  ) ! Size of Krylov Space'  >>SIZE
  fi
  cat SIZE | grep -i 'cv_delt' >/dev/null
  if [ $? -ne 0 ]; then
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      real cv_delt' >> SIZE
      echo '      parameter(cv_delt = 0.1 ) ! linear convergence factor ' >>SIZE
  fi
  cat SIZE | grep -i 'cv_sigs' >/dev/null
  if [ $? -ne 0 ]; then
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      real cv_sigs' >> SIZE
      echo '      parameter(cv_sigs = 0.01) ! Jacobian DQ perturbation scaling factor'  >>SIZE
  fi
  fi
  cat SIZE | grep -i 'toteq' >/dev/null
  if [ $? -ne 0 ]; then
  if [ "$IFCMT" == "true" ]; then
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      integer toteq' >> SIZE
      echo '      parameter(toteq = 5  ) ! Number of conserved variables '  >>SIZE
      echo 'c IFCMT=TRUE  then toteq=5'  >>SIZE
  else
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      integer toteq' >> SIZE
      echo '      parameter(toteq = 1  ) ! Number of conserved variables '  >>SIZE
      echo 'c IFCMT=FALSE  then toteq=1'  >>SIZE
  fi
  fi
  cat SIZE | grep -i 'lelcmt' >/dev/null
  if [ $? -ne 0 ]; then
  if [ "$IFCMT" == "true" ]; then
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      integer lelcmt' >> SIZE
      echo '      parameter(lelcmt = lelt  ) ! # of cmt elements '  >>SIZE
      echo 'c IFCMT=TRUE  then lelcmt=lelt'  >>SIZE
  else
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      integer lelcmt' >> SIZE
      echo '      parameter(lelcmt = 1  ) ! # of cmt elements '  >>SIZE
      echo 'c IFCMT=FALSE  then lelcmt=1'  >>SIZE
  fi
  fi
else
  echo "FATAL ERROR: Cannot find SIZE"
  exit 1
fi

# tweak .usr (actually, .f ) file
rm -f $CASENAME.f
cp -p $CASENAME.usr $CASENAME.f
cat $CASENAME.f | grep -i "subroutine.*usrsetvert" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
_ACEOF
fi

if [ "$IFCMT" == "true" -o "$IFCMT" == "yes" ]; then
cat $CASENAME.f | grep -i "subroutine.*cmt_switch" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine cmt_switch ! to set IFCMT logical flag
      include 'SIZE'
      include 'INPUT'
      IFCMT=.true.
      return
      end
_ACEOF
fi
else
cat $CASENAME.f | grep -i "subroutine.*cmt_switch" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine cmt_switch ! to set IFCMT logical flag
      include 'SIZE'
      include 'INPUT'
      IFCMT=.false.
      return
      end
_ACEOF
fi
fi

cat $CASENAME.f | grep -i "subroutine.*usrflt" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine usrflt(rmult) ! user defined filter
      include 'SIZE'
      real rmult(lx1)
      call rone(rmult,lx1)
      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*userflux" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine userflux ! user defined flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real fluxout(lx1*lz1)
      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*userEOS" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine userEOS ! user defined EOS 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*userq2" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
      subroutine userq2(q) ! set source term w/ MOAB

#include "NEKMOAB"

      include 'GEOM'
      include 'TSTEP'

      common /usrparl/ ifnekalone,mbpcid
      logical          ifnekalone   ! true for compatibility when p123>0
      common /moabqv/   qvola(lx1,ly1,lz1,lelt)     ! vol. heat s. for MOAB
      integer          mbpcid

      real q(lx1*ly1*lz1,lelt)

      integer e, i, ime, ntot, nxyz
      real x, y, r, init_power_density

      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelt

      if (istep .eq. 0 .or. (.not.ifcoup .and. .not.ifvcoup) ) then

         init_power_density = param(121)
         call rzero (q,ntot)
         if (ifnekalone) then
           do e = 1,nelt
             do i = 1,nxyz
               x = xm1(i,1,1,e)
               y = ym1(i,1,1,e)
               r = (x-16.d0)**2 + (y-7.4019221593017583d0)**2  ! heat source region
               r = sqrt(r)

               if (2.0.le.r .and. r.le.4.0) q(i,e) = init_power_density   ! w/ clad
             enddo
           enddo
         else
           do e = 1,nelt
             do i = 1,nxyz

c .rea:         4 2 3 5 6       fluid set #4, solid set #2 3 5 6
c               1 2 3 5 6       ! imatie: fluid fuel clad crod can

               ime = imatie(e)

               if (ime.eq.2) q(i,e) = init_power_density  ! heat source in the solid -- w/ clad
             enddo
           enddo
         endif

      else

         if (ifcoup) then
            call nekMOAB_import_vars(powTag, 0, q)
         elseif (ifvcoup) then
            call nekMOAB_import_vars(vpowTag, 1, q)
         elseif (istep.lt.10) then
           write(6,*) "NEKMOAB:userq2: Not importing power density."
           write(6,*) "NEKMOAB:userq2 *** Enable IFCOUP or IFVCOUP ***"
         endif

      endif

      call compute_total_power(q)

      call copy(qvola,q,ntot)   ! for output

      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*compute_total_power" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
      subroutine compute_total_power(q)
      include 'SIZE'
      include 'TOTAL'

      integer ntot
      real q(lx1*ly1*lz1,lelt), bmq(lx1*ly1*lz1*lelt)

      ntot = nx1*ny1*nz1*nelt          ! Total power calculation
      call col3(bmq,q,bm1,ntot)
      qtot = glsum(bmq,ntot)
      if (nid.eq.0)  write(6,1) istep,time,qtot
    1 format (i8,1p2e15.7,' W  Qtot')

      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*usr_get_current_time" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
      subroutine usr_get_current_time(current_time)
      include 'SIZE'
      include 'TSTEP'
      real current_time
      
      integer ispert, glob_tstep
      common /usrperturbation/  ispert, glob_tstep  ! perturbation of temperature and density BC
      
      if (ispert .eq. 0) then
        current_time = 0.0
      else
        current_time = TIME
      endif
      
      print *, "[COUPE:NEK]...Getting current final time from NEK ",
     $                        current_time

      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*nekMOAB_reimport" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
      subroutine nekMOAB_reimport
C      implicit none
C  include "NEKMOAB"
C      include 'PARALLEL'
C      include 'GEOM'
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'ZPER'
      include 'OPCTR'

      integer ierr, igeom

      call mapelpr                        ! create gllel mapping 
      call moab_geometry(xm1,ym1,zm1)     ! fill xm1,ym1,zm1
      call xml2xc                         ! fill xc,yc,zc

C      call nekMOAB_BC             ! read MOAB BCs 

C      call nekMOAB_compute_diagnostics

C     generate geometry (called after usrdat in case something changed)
      igeom=2
      call gengeom (igeom)

      if (ifmvbd) call setup_mesh_dssum ! Set up dssum for mesh (needs geom)

      call geom_reset(igeom)    ! recompute Jacobians, etc.

      call vrdsmsh          ! verify mesh topology

      return
      end  
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*usr_reset_pseudo_time" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
      subroutine usr_reset_pseudo_time
      include 'SIZE'
      include 'TSTEP'

      ISTEP = 0
      LASTEP = 0
C      NSTEPS = 10

      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*usr_set_pcid" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
      subroutine usr_set_pcid(ipid)

      common /usrparl/ ifnekalone,mbpcid
      logical          ifnekalone   ! true for compatibility when p123>0
      integer          mbpcid,ipid

      mbpcid=ipid

      return
      end
_ACEOF
fi

cat $CASENAME.f | grep -i "subroutine.*usr_moab_output" >/dev/null
if [ $? -ne 0 ]; then
cat >> $CASENAME.f << _ACEOF
c
c automatically added by SHARP makenek
#ifdef MOAB
      subroutine usr_moab_output()

      implicit none
#include "NEKMOAB"
      include 'TSTEP'
      include 'PARALLEL'

      integer ierr
      CHARACTER*20 NAME, STEPSTR
      CHARACTER*132 wroptions, seroptions, paroptions
      CHARACTER*1  NUMRL(0:9)
      DATA NUMRL          /'0','1','2','3','4','5','6','7','8','9'/
      data seroptions /''/
      data paroptions /' moab:PARALLEL=WRITE_PART moab:PCOMM=1 '/
c      data paroptions /' moab:PARALLEL=WRITE_PART moab:PCOMM=0 '/
      
      integer ispert, glob_tstep
      common /usrperturbation/  ispert, glob_tstep  ! perturbation of temperature and density BC
      
      common /usrparl/ ifnekalone,mbpcid
      logical          ifnekalone   ! true for compatibility when p123>0
      integer          mbpcid

      if (np .gt. 1) then
        mbpcid=1
        if (mbpcid .ge. 0) then
          WRITE(wroptions,'(A,A,I1)')' moab:PARALLEL=WRITE_PART ',
     $        'moab:PCOMM=',mbpcid
        else
          wroptions=paroptions
        endif
        wroptions=TRIM(wroptions)
      else
        wroptions=TRIM(seroptions)
      endif

c      WRITE(STEPSTR, '(I5.5)') glob_tstep
      WRITE(STEPSTR,'(A,I5.5,A)') 'sahex',glob_tstep,'.h5m'
      NAME=TRIM(STEPSTR)

      if (ispert .eq. 0) glob_tstep = glob_tstep + 1

c copy data from nek to MOAB
      call nekMOAB_export_vars

      ierr = 0
      if(istep .eq. 0 .or. iostep .eq. 0 .or.
     $     mod(istep,iostep) .eq. 0) then
         if (np .gt. 1) then
c            call iMeshP_saveAll(%VAL(imeshh),%VAL(hPartn),%VAL(fileset),
c     $           name,wroptions,ierr)
         else
            call iMesh_save(%VAL(imeshh),%VAL(fileset),name,wroptions,
     $           ierr)
         endif
c no check for error, just print something and continue
         if (ierr .ne. 0 .and. nid .eq. 0) then
           write(6,*) 'Step:',istep,' Writing to ',
     $                   TRIM(name),' = failed. ierr=',ierr
         else if (nid .eq. 0) then
			     write(6,*) 'Step:',istep,' Writing to ',
     $                   TRIM(name),' = success'
         endif
      endif

      return
      end

#endif
_ACEOF
fi

if [ $NOCOMPILE -eq 1 ]; then
  exit 0
fi 
