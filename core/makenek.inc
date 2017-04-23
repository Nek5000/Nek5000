# This include file is used by the makenek script
# to automatically create a makefile for Nek5000 
# ------------------------------------------------

echo "makenek - automatic build tool for Nek5000"

if [ "$PPLIST" == "?" ]; then
  echo "available pre-processor symbols:" 
  echo "  NOMPIIO     deactivate MPI-IO support"
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
  echo "usage: makenek [.usr filename | clean | realclean] < -nocompile >"
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

if [ $1 = "realclean" ]; then
  make -C $SOURCE_ROOT_JL clean
  make clean
  rm -rf makefile
  exit 0
fi

NOCOMPILE=0
if [ "$2" == "-nocompile" ]; then
  NOCOMPILE=1
fi 

CASENAME=$1
CASEDIR=`pwd`

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
   mv -v SIZEu SIZE
fi

if [ ! -f SIZE ]; then
   echo "FATAL ERROR: Cannot find SIZE!"
   exit 1
fi

if [ ! -f $SOURCE_ROOT/core/makefile.template ]; then
  echo "FATAL ERROR: Cannot find $SOURCE_ROOT/core/makefile.template!"
  exit 1
fi

# test if F77 compiler exists
if ! which `echo $F77 | awk '{print $1}'` 1>/dev/null ; then
  echo "FATAL ERROR: Cannot find `echo $F77 | awk '{print $1}'`!"
  exit 1
fi

# test if C compiler exists
if ! which `echo $CC | awk '{print $1}'` 1>/dev/null ; then 
  echo "FATAL ERROR: Cannot find `echo $CC | awk '{print $1}'`!"
  exit 1
fi

# initial clean-up
rm -f nek5000 2>/dev/null
rm -f ./obj/subuser.o        2>/dev/null

# Check if the compiler adds an underscore to external functions
UNDERSCORE=0
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call test
      end
_ACEOF
$F77 -c test_underscore.f 2>&1 >/dev/null 
if nm test_underscore.o | grep -q test_ ; then
  UNDERSCORE=1
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
  *pgf*)        FFLAGS+=" -r8 -Mpreprocess"
               ;;
  *gfortran*)   FFLAGS+=" -fdefault-real-8 -fdefault-double-8 -cpp"
               ;;
  *ifort*)      FFLAGS+=" -r8 -fpconstant -fpp"
               ;;
  *pathf*)      FFLAGS+=" -r8 -cpp -fno-second-underscore"
               ;;
  *ftn*)        FFLAGS+=" -s default64 -eF"
               ;;
  *xlf*)       FFLAGS+=" -qrealsize=8 -qdpc=e -qsuffix=cpp=f"
               PPPO="-WF,"
               F77="${F77} -qsuppress=cmpmsg"
               ;;
  *)  echo "ERROR: Unable to detect compiler!"
      echo "       Please edit the makefile and specify FFLAGS and CFLAGS"
      echo ""
      NOCOMPILE=1
      read;;
esac

# set preprocessor symbols
if [ $MPI -eq 0 ]; then
  PPLIST="${PPLIST} NOMPIIO"
else
  MPI=1
  PPLIST="${PPLIST} MPI"
fi

if [ $UNDERSCORE -ne 0 ]; then
  PPLIST="${PPLIST} UNDERSCORE"
fi 

if [ $PROFILING -ne 0 ]; then
  if [ $MPI -ne 0 ]; then
    PPLIST="${PPLIST} TIMER"
  else
    PPLIST="${PPLIST} TIMER"
  fi
fi

PPLIST="${PPLIST}"

# add user defined preprocessor symbols
CMT=0
NEKNEK=0
for i in $PPLIST
do
   if [ "$i" == "NEKNEK" ]; then
     NEKNEK=1
   fi
   if [ "$i" == "CMTNEK" ]; then
     CMT=1
   fi
done

if echo $PPLIST | grep -q 'CVODE' ; then 
   CVODE=1 
else
   CVODE=0
fi

MXM_USER="mxm_std.o"
if echo $PPLIST | grep -q 'BGQ' ; then 
   MXM_USER+=" mxm_bgq.o" 
fi

BLAS="blas.o dsygv.o"
if echo $PPLIST | grep -q 'VENDOR_BLAS' ; then 
   BLAS=" "
fi

# tweak makefile template 
echo "generating makefile ..."
rm -rf makefile 2>/dev/null

sed -e "s:^FC[ ]*=.*:FC\:=$F77:" \
-e "s:^CC[ ]*=.*:CC\:=$CC:" \
-e "s:^FFLAGS[ ]*=.*:FFLAGS\:=$FFLAGS:" \
-e "s:^CFLAGS[ ]*=.*:CFLAGS\:=$CFLAGS:" \
-e "s/^PPPO[ ]*=.*/PPPO=$PPPO/" \
-e "s/^PPS[ ]*=.*/PPS=$PPLIST/" \
-e "s:^MXM[ ]*=.*:MXM=$MXM_USER:" \
-e "s:^BLAS[ ]*=.*:BLAS=$BLAS:" \
-e "s/^MPI[ ]*=.*/MPI:=$MPI/" \
-e "s/^CMT[ ]*=.*/CMT:=$CMT/" \
-e "s/^NEKNEK[ ]*=.*/NEKNEK:=$NEKNEK/" \
-e "s/^VISIT[ ]*=.*/VISIT:=${VISIT}/" \
-e "s:^USR[ ]*=.*:USR\:=$USR:" \
-e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS\:=$USR_LFLAGS:" \
-e "s:^S[ ]*=.*:S\:=${SOURCE_ROOT}/core:" \
-e "s:^SCMT[ ]*=.*:SCMT\:=${SOURCE_ROOT}/core/cmt:" \
$SOURCE_ROOT/core/makefile.template >.makefile

if [ $VISIT -ne 0 ]; then
  echo ""
  echo "Activate VisIt Sim Mode"

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

  sed -e "s/^VISIT[ ]*=.*/VISIT:=${VISIT}/" \
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
  if [ $CMT -ne 0 ]; then
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


if [ $CVODE -ne 0 ]; then

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

if [ $CMT -ne 0 ]; then
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

export NEK_WARN 
export MPI
export FFLAGS
export CFLAGS
export UNDERSCORE

if [ $NOCOMPILE -ne 0 ]; then
  exit 0
fi 
