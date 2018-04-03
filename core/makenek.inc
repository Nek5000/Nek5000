function make_3rd_party() {
  if [ $NOCOMPILE -ne 0 ]; then
    return 0
  fi
  echo "build 3rd-party software ..."

  cd $SOURCE_ROOT/3rd_party/gslib
  ./install
  cd $CASEDIR 

  if [ $CVODE -ne 0 ]; then
     cd $SOURCE_ROOT/3rd_party/cvode
     ./install
  fi
  cd $CASEDIR
}

function clean_3rd_party() {
  echo "cleaning 3rd-party software ..."
  cd $SOURCE_ROOT/3rd_party/gslib; ./install clean 2>/dev/null
  cd $CASEDIR
  cd $SOURCE_ROOT/3rd_party/cvode; ./install clean 2>/dev/null
  cd $CASEDIR
}

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
  echo "  CMTNEK      activate DG compressible-flow solver (experimental)"
  exit 1
fi

# first do some checks ...
if [ "$1" == "-h" ]; then
  echo ""
  echo "usage: makenek [.usr filename | clean] < -nocompile >"
  echo ""
  exit 1
fi

CASEDIR=`pwd`
CASENAME=$1
if [ $# -eq 0 ]; then
   NCASEFILES=`find *.usr | wc -l`
   if [ $NCASEFILES != "1" ]; then
     echo "FATAL ERROR: Found multiple .usr files!"
     echo "Please use makenek [.usr filename]"
     exit 1
   else
     CASENAME=`find *.usr | sed 's/\.usr$//g'`
   fi
fi

if [ "$mver" != "17.0.4" ]; then
  echo "FATAL ERROR: Your makenek is incompatible!"
  echo "Please upgrade to the version shipped with your Nek source."
  exit 1
fi

if [ ! "$SOURCE_ROOT_GSLIB" ]; then
  SOURCE_ROOT_GSLIB="$SOURCE_ROOT/3rd_party/gslib/lib"
  export SOURCE_ROOT_GSLIB
fi

if [[ $1 = "clean" ]]; then
  while true; do
    read -p "do you want to clean all 3rd party software? [N]" yn
    case $yn in
      [Yy]* )
        clean_3rd_party; 
        make clean;
        break ;;
      * ) 
        make clean;
        break ;;
    esac
  done
  rm -rf compiler.out
  rm -rf makefile
  rm -rf .state
  exit 0
fi

NOCOMPILE=0
if [[ "$2" == "-nocompile" ]]; then
  NOCOMPILE=1
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

# test if FC compiler exists
if ! which `echo $FC | awk '{print $1}'` 1>/dev/null ; then
  echo "FATAL ERROR: Cannot find `echo $FC | awk '{print $1}'`!"
  exit 1
fi

# test if C compiler exists
if ! which `echo $CC | awk '{print $1}'` 1>/dev/null ; then 
  echo "FATAL ERROR: Cannot find `echo $CC | awk '{print $1}'`!"
  exit 1
fi

# trying to figure which compiler the wrapper is using 
FCok=0

FCcomp_=`$FC -showme 2>/dev/null | head -1 2>/dev/null 1>.tmp || true`
FCcomp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}' || true`
if [ -f "`which $FCcomp 2>/dev/null`" ]; then
  FCok=1
fi

if [ $FCok -eq 0 ]; then
  FCcomp_=`$FC -show 2>/dev/null | head -1 2>/dev/null 1>.tmp || true`
  FCcomp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}' || true`
  if [ -f "`which $FCcomp 2>/dev/null`" ]; then
    FCok=1
  fi
fi

if [ $FCok -eq 0 ]; then
  FCcomp_=`$FC -craype-verbose 2>/dev/null 1>.tmp || true`
  FCcomp=`cat .tmp | awk '{print $1}' | awk -F/ '{print $NF}' || true`
  if [ -f "`which $FCcomp 2>/dev/null`" ]; then
    FCok=1
  fi
fi

if [ $FCok -eq 0 ]; then
  FCcomp=`echo $FC | awk '{print $1}'`
  if [ -f "`which $FCcomp 2>/dev/null`" ]; then
    FCok=1
  fi
fi

\rm -f .tmp
if [ $FCok -eq 0 ]; then
  FCcomp="unknown"
fi

# assign FC compiler specific flags
case $FCcomp in
  *pgf*)        FFLAGS+=" -r8 -Mpreprocess"
               ;;
  *gfortran*)   FFLAGS+=" -fdefault-real-8 -fdefault-double-8 -cpp"
               ;;
  *ifort*)      FFLAGS+=" -r8 -fpconstant -fpp"
               ;;
  *xlf*)        FFLAGS+=" -qfixed -qrealsize=8 -qdpc=e -qsuffix=cpp=f -qsuppress=cmpmsg"
                PPPO="-WF,"
               ;;
  *)  echo "ERROR: Cannot find a supported compiler!"
      echo ""
      NOCOMPILE=1
      read;;
esac

# initial clean-up
rm -f $CASENAME.f nek5000 2>/dev/null
rm -f ./obj/subuser.o 2>/dev/null

# Check if the compiler adds an underscore to external functions
UNDERSCORE=0
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call test
      end
_ACEOF
$FCcomp -c test_underscore.f 2>&1 >/dev/null 
if nm test_underscore.o | grep -q test_ ; then
  UNDERSCORE=1
fi
\rm test_underscore.* 2>/dev/null

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
MPIIO=1
for i in $PPLIST
do
   if [ "$i" == "CMTNEK" ]; then
     CMT=1
   fi
   if [ "$i" == "NOMPIIO" ]; then
     MPIIO=0
   fi
done

if echo $PPLIST | grep -q 'CVODE' ; then 
   CVODE=1
   if [ ! "$LIBPATH_CVODE" ]; then
     LIBPATH_CVODE=$SOURCE_ROOT/3rd_party/cvode/lib 
   fi
   if [ $MPI -eq 1 ]; then
     USR_LFLAGS+=" -L$LIBPATH_CVODE -lsundials_fcvode -lsundials_cvode"
     USR_LFLAGS+=" -lsundials_fnvecparallel -lsundials_nvecparallel" 
   else
     USR_LFLAGS+=" -L$LIBPATH_CVODE -lsundials_fcvode -lsundials_cvode"
     USR_LFLAGS+=" -lsundials_fnvecserial -lsundials_nvecserial" 
   fi
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

# gslib build options
GSLIB_OPT+=" MPI=$MPI MPIIO=$MPIIO"
GSLIB_OPT+=" ADDUS=$UNDERSCORE USREXIT=1 BLAS=2"
GSLIB_OPT+=" PREFIX=gslib_ FPREFIX=fgslib_"
USR_LFLAGS+=" -L$SOURCE_ROOT_GSLIB -lgs"

# tweak makefile template 
echo "generating makefile ..."
rm -rf makefile 2>/dev/null

sed -e "s:^FC[ ]*=.*:FC\:=$FC:" \
-e "s:^CC[ ]*=.*:CC\:=$CC:" \
-e "s:^FFLAGS[ ]*=.*:FFLAGS\:=$FFLAGS:" \
-e "s:^CFLAGS[ ]*=.*:CFLAGS\:=$CFLAGS:" \
-e "s/^PPPO[ ]*=.*/PPPO=$PPPO/" \
-e "s/^PPS[ ]*=.*/PPS=$PPLIST/" \
-e "s:^MXM[ ]*=.*:MXM=$MXM_USER:" \
-e "s:^BLAS[ ]*=.*:BLAS=$BLAS:" \
-e "s/^MPI[ ]*=.*/MPI:=$MPI/" \
-e "s/^CVODE[ ]*=.*/CVODE:=$CVODE/" \
-e "s/^CMT[ ]*=.*/CMT:=$CMT/" \
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

# tweak leagacy SIZE file
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
  if ! cat SIZE | grep -qi 'nsessmax' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (nsessmax=1)  ! max sessions to NEKNEK' >> SIZE
  fi
  if ! cat SIZE | grep -qi 'nmaxl_nn' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (nmaxl_nn=' >> SIZE
     echo '     $          min(1+(nsessmax-1)*2*ldim*lxz*lelt,2*ldim*lxz*lelt))' >>SIZE
  fi
  if ! cat SIZE | grep -qi 'nfldmax_nn' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      parameter (nfldmax_nn=' >> SIZE
     echo '     $          min(1+(nsessmax-1)*(ldim+1+ldimt),ldim+1+ldimt))' >>SIZE
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

# tweak all SIZE file versions
if ! cat SIZE | grep -qi 'ldimt_proj' ; then
     echo >>SIZE
     echo 'c automatically added by makenek' >>SIZE
     echo '      integer ldimt_proj' >>SIZE
     echo '      parameter(ldimt_proj = 1) ! max auxiliary fields residual projection ' >>SIZE
fi

export FC
export CC
export NEK_WARN 
export MPI
export MPIO
export FFLAGS
export CFLAGS
export UNDERSCORE
export GSLIB_OPT


if [ $NOCOMPILE -ne 0 ]; then
  exit 0
fi

declare -a varlist=(
"SOURCE_ROOT"
"FC"
"FFLAGS"
"CC"
"CFLAGS"
"PPLIST"
)
if [ -f .state ]; then
  for i in "${varlist[@]}"
  do
    if ! cat .state | grep -q "$i=${!i}" ; then
      echo "configuration change detected, forcing full rebuild"
      rm -rf .state
      while true; do
        read -p "do you want to rebuild all 3rd party software? [N]" yn
        case $yn in
          [Yy]* )
            clean_3rd_party;
            break ;;
          * )
           break ;;
        esac
      done  
      make clean;
      break
    fi
  done ;
fi

rm -rf .state
for i in "${varlist[@]}"
do
  echo "$i=${!i}" >> .state
done ;
