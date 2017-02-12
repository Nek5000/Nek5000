# This include file is used by the makenek script
# to automatically create a makefile for Nek5000 
# ------------------------------------------------

echo "makenek - automatic build tool for Nek5000"

if [ "$PPLIST" == "?" ]; then
  echo "available pre-processor symbols:" 
  echo "  AMG         use AMG instead of XXT solver"
  echo "  MPIIO       use MPI-IO for I/O"
  echo "  BGQ         use BGQ optimized mxm"
  echo "  XSMM        use libxsmm for mxm"
  echo "  CVODE       compile with CVODE support for scalars"
  echo "  VENDOR_BLAS use VENDOR BLAS/LAPACK"
  echo "  EXTBAR      add underscore to exit call (for BGQ)"
  echo "  NEKNEK      activate overlapping mesh solver (experimental)"
  echo "  CMTNEK      activate DG compressible-flow solver (experimental)"
  exit 1
fi

# first do some checks ...
if [ $# -eq 0 ]; then
  echo ""
  echo "usage: makenek [.usr filename | clean] < -nocompile >"
  echo ""
  exit 1
fi

if [ "$mver" != "17.0.0" ]; then
  echo "ERROR: Your makenek is incompatible!"
  echo "Please upgrade to the version shipped with your Nek source."
  exit 1
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
APATH_SRC=`cd $SOURCE_ROOT_CORE; pwd`
SOURCE_ROOT_CORE=$APATH_SRC
APATH_SRC=`cd $SOURCE_ROOT_JL; pwd`
SOURCE_ROOT_JL=$APATH_SRC
APATH_SRC=`cd $SOURCE_ROOT_CMT; pwd`
SOURCE_ROOT_CMT=$APATH_SRC

IFNEKNEK=false

# do some basic checks
if [ "$CASEDIR" == "$SOURCE_ROOT_CORE" ]; then
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

if [ ! -f $SOURCE_ROOT_CORE/makefile.template ]; then
  echo "FATAL ERROR: Cannot find $SOURCE_ROOT_CORE/makefile.template!"
  exit 1
fi

# test if F77 compiler exists
if ! which -s `echo $F77 | awk '{print $1}'` ; then
  echo "FATAL ERROR: Cannot find `echo $F77 | awk '{print $1}'`!"
  exit 1
fi

# test if C compiler exists
if ! which -s `echo $CC | awk '{print $1}'` ; then 
  echo "FATAL ERROR: Cannot find `echo $CC | awk '{print $1}'`!"
  exit 1
fi

# initial clean-up
rm -f nek5000 2>/dev/null
rm -f ./obj/subuser.o        2>/dev/null

# Check if the compiler adds an underscore to external functions
UNDERSCORE=false
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call test
      end
_ACEOF
$F77 -c test_underscore.f 2>&1 >/dev/null 
if nm test_underscore.o | grep -q test_ ; then
  UNDERSCORE=true
fi
\rm test_underscore.* 2>/dev/null


# trying to figure which compiler the wrapper is using 
F77ok=0

F77comp_=`$F77 -showme 2>/dev/null 1>.tmp || true`
F77comp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}' || true`
if [ -f "`which $F77comp 2>/dev/null`" ]; then
  F77ok=1
fi

if [ $F77ok -eq 0 ]; then
  F77comp_=`$F77 -show 2>/dev/null 1>.tmp || true`
  F77comp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}' || true`
  if [ -f "`which $F77comp 2>/dev/null`" ]; then
    F77ok=1
  fi
fi

if [ $F77ok -eq 0 ]; then
  F77comp_=`$F77 -craype-verbose 2>/dev/null 1>.tmp || true`
  F77comp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}' || true`
  if [ -f "`which $F77comp 2>/dev/null`" ]; then
    F77ok=1
  fi
fi

if [ $F77ok -eq 0 ]; then
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
  *pgf*)        FFLAGS="-r8 -Mpreprocess"
               ;;
  *gfortran*)   FFLAGS="-fdefault-real-8 -fdefault-double-8 -x f77-cpp-input"
               ;;
  *ifort*)      FFLAGS="-r8 -fpconstant -fpp"
               ;;
  *pathf*)      FFLAGS="-r8 -cpp -fno-second-underscore"
               ;;
  *ftn*)        FFLAGS="-s default64 -eF"
               ;;
  *xlf*)       FFLAGS="-qrealsize=8 -qdpc=e -qsuffix=cpp=f"
               PPPO="-WF,"
               F77="${F77} -qsuppress=cmpmsg"
               ;;
  *)  echo "ERROR: Unable to detect compiler!"
      echo "       Please edit the makefile and specify FFLAGS and CFLAGS"
      echo ""
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

if [ "$UNDERSCORE" == "true" ]; then
  PPLIST="${PPLIST} UNDERSCORE"
fi 

if [ "$IFPROFILING" == "true" -o "$IFPROFILING" == "yes" ]; then
  if [ "$IFMPI" == "true" -o "IFMPI" == "yes" ]; then
    PPLIST="${PPLIST} TIMER MPITIMER"
  else
    PPLIST="${PPLIST} TIMER"
  fi
fi

PPLIST="${PPLIST} GLOBAL_LONG_LONG"

# add user defined preprocessor symbols
IFCMT=false
for i in $PPLIST
do
   if [ "$i" == "NEKNEK" ]; then
     IFNEKNEK=true
   fi
   if [ "$i" == "CMTNEK" ]; then
     IFCMT=true
   fi
done
export IFCMT

if echo $PPLIST | grep -q 'CVODE' ; then 
   IFCVODE=true 
else
   IFCVODE=false
fi

# AMG setup
if echo $PPLIST | grep -q 'AMG' ; then 
   IFAMG=true 
else
   IFAMG=false
fi

MXM_USER="mxm_std.o"
if echo $PPLIST | grep -q 'BGQ' ; then 
   MXM_USER+=" mxm_bgq.o" 
fi

BLAS="blas.o dsygv.o ssygv.o"
if echo $PPLIST | grep -q 'VENDOR_BLAS' ; then 
   BLAS=" "
fi

if [ "$USR_LIB" != "" ]; then
    USR_LFLAGS="${USR_LFLAGS} ${USR_LIB}"
fi

# tweak makefile template 
echo "generating makefile ..."
rm -rf makefile 2>/dev/null

sed -e "s:^FC[ ]*=.*:FC\:=$F77:" \
-e "s:^CC[ ]*=.*:CC\:=$CC:" \
-e "s:^G[ ]*=.*:G\:=$G:" \
-e "s/^FFLAGS[ ]*=.*/FFLAGS:=$FFLAGS/" \
-e "s/^CFLAGS[ ]*=.*/CFLAGS:=$CFLAGS/" \
-e "s/^PPPO[ ]*=.*/PPPO=$PPPO/" \
-e "s/^PPS[ ]*=.*/PPS=$PPLIST/" \
-e "s:^MXM[ ]*=.*:MXM=$MXM_USER:" \
-e "s:^BLAS[ ]*=.*:BLAS=$BLAS:" \
-e "s/^IFAMG[ ]*=.*/IFAMG:=$IFAMG/" \
-e "s/^IFMPI[ ]*=.*/IFMPI:=$IFMPI/" \
-e "s/^IFCMT[ ]*=.*/IFCMT:=$IFCMT/" \
-e "s/^IFNEKNEK[ ]*=.*/IFNEKNEK:=$IFNEKNEK/" \
-e "s:^USR[ ]*=.*:USR\:=$USR:" \
-e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS\:=$USR_LFLAGS:" \
-e "s:^S[ ]*=.*:S\:=${SOURCE_ROOT_CORE}:" \
-e "s:^S2[ ]*=.*:S2\:=${SOURCE_ROOT_CMT}:" \
-e "s:^J[ ]*=.*:J\:=${SOURCE_ROOT_JL}:" $SOURCE_ROOT_CORE/makefile.template >.makefile

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
if ! cat SIZE | grep -q 'SIZE.inc' ; then
  if ! cat SIZE | grep -qi 'optlevel' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      integer optlevel,loglevel' >>SIZE
     echo '      common /lolevels/ optlevel,loglevel' >>SIZE
  fi
  if ! cat SIZE | grep -qi 'lxo' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)' >>SIZE
  fi
  if ! cat SIZE | grep -qi 'ax1' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      integer ax1,ay1,az1,ax2,ay2,az2' >> SIZE
     echo '      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! running averages' >> SIZE
  fi
  if ! cat SIZE | grep -qi 'lys=lxs' ; then
     cat SIZE | grep -iv lxs > SIZE.x; mv SIZE.x SIZE  # Clean existing SIZE file of old version
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Preconditioner' >> SIZE

  fi
  if ! cat SIZE | grep -qi 'lcvx1' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      integer lcvx1,lcvy1,lcvz1,lcvelt' >> SIZE
     echo '      parameter (lcvx1=1,lcvy1=1,lcvz1=1,lcvelt=1) ! cvode arrays' >> SIZE
  fi
  if ! cat SIZE | grep -qi 'lfdm' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (lfdm=0)  ! == 1 for fast diagonalization method' >> SIZE
  fi
  if ! cat SIZE | grep -qi 'nio' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      common/IOFLAG/nio  ! for logfile verbosity control' >> SIZE
  fi

  if ! cat SIZE | grep -qi 'toteq' ; then
  if [ "$IFCMT" == "true" ]; then
      echo >>SIZE
      echo 'c automatically added by makenek' >>SIZE
      echo '      integer toteq' >> SIZE
      echo '      parameter(toteq = 5  ) ! Number of conserved variables '  >>SIZE
  fi
  fi
fi

# tweak .usr (actually, .f ) file
rm -f $CASENAME.f
cp -pv $CASENAME.usr $CASENAME.f

if ! cat $CASENAME.f | grep -qi "subroutine.*usrsetvert" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end
_ACEOF
fi

if ! cat $CASENAME.f | grep -qi "subroutine.*userqtl" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
_ACEOF
fi

if [ "$IFCVODE" == "true" -o "$IFCVODE" == "yes" ]; then

if ! cat $CASENAME.f | grep -qi "^#include.*cvode_aux.*\.h" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
#include "cvode_aux.h"
_ACEOF
fi

if ! cat $CASENAME.f | grep -qi "^#include.*cvode_preco.*\.h" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
#include "cvode_preco_dummy.h"
_ACEOF
fi

fi

if [ "$IFCMT" == "true" -o "$IFCMT" == "yes" ]; then
if ! cat $CASENAME.f | grep -qi "subroutine.*cmt_usrflt" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
      subroutine cmt_usrflt(rmult) ! user defined filter
      include 'SIZE'
      real rmult(lx1)
      call rone(rmult,lx1)
      return
      end
_ACEOF
fi

if ! cat $CASENAME.f | grep -qi "subroutine.*cmt_userflux" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
      subroutine cmt_userflux ! user defined flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTDATA'
      real fluxout(lx1*lz1)
      return
      end
_ACEOF
fi

if ! cat $CASENAME.f | grep -qi "subroutine.*cmt_userEOS" ; then
cat >> $CASENAME.f << _ACEOF

c automatically added by makenek
      subroutine cmt_userEOS ! user defined EOS 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTDATA'

      return
      end
_ACEOF
fi
if [ -e "cmtparticles.usrp" ]
then
echo 'Particles found CMT'
else
cat >> $CASENAME.f << _ACEOF
c
c automatically added by makenek
      subroutine usr_particles_init ! used for particles
      return
      end
c
c automatically added by makenek
      subroutine usr_particles_solver ! used for particles
      return
      end
c
c automatically added by makenek
      subroutine usr_particles_io(istep) ! used for particles
      integer istep
      return
      end

_ACEOF
fi
fi


if [ $NOCOMPILE -eq 1 ]; then
  exit 0
fi 
