trap onexit EXIT

is_llvm_flang() {
    local cmd="$1"

    if ! command -v "$cmd" &> /dev/null; then
        return 1
    fi

    if "$cmd" -fc1 -help >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

function onexit() {
  retval=$?
  if [ $retval -gt 0 ] && [ $retval -lt 11 ]; then
     printf "\nAn ERROR occured, see build.log for details\n"
  fi
}

function build_3rdparty() {
  printf "building 3rd-party dependencies ... "
  set -o pipefail

  cd $SOURCE_ROOT_BLAS
  FFLAGS="-O0 $FF77 $FFLAGS" ./install 2>&1 | tee -a $CASEDIR/build.log >/dev/null
  cd $CASEDIR 

  cd $SOURCE_ROOT_GSLIB
  CFLAGS="-O2 $CFLAGS" ./install 2>&1 | tee -a $CASEDIR/build.log >/dev/null
  cd $CASEDIR 

  if [ $PARRSB -ne 0 ]; then
     cd $SOURCE_ROOT_PARRSB
     CFLAGS="-O2 $CFLAGS" ./install 2>&1 | tee -a $CASEDIR/build.log >/dev/null
  fi
  cd $CASEDIR

  if [ $CVODE -ne 0 ]; then
     cd $SOURCE_ROOT_CVODE
     CFLAGS="-O2 $CFLAGS" ./install 2>&1 | tee -a $CASEDIR/build.log >/dev/null
  fi
  cd $CASEDIR

  if [ $HYPRE -ne 0 ]; then
     cd $SOURCE_ROOT_HYPRE
     CFLAGS="-O2 $CFLAGS" ./install 2>&1 | tee -a $CASEDIR/build.log >/dev/null
  fi
  cd $CASEDIR

  printf "done\n"
  set +o pipefail
}

function clean_3rdparty() {
  if [[ $1 = "clean" ]]; then
    act="clean"
  elif [[ $1 = "distclean" ]]; then
    act="distclean"
  else
    echo "clean_3rdparty: Invalid clean mode $1"
    exit 1
  fi

  cd $SOURCE_ROOT_BLAS; ./install $act 2>/dev/null
  cd $CASEDIR
  cd $SOURCE_ROOT_GSLIB; ./install $act 2>/dev/null
  cd $CASEDIR
  cd $SOURCE_ROOT_PARRSB; ./install $act 2>/dev/null
  cd $CASEDIR
  cd $SOURCE_ROOT_CVODE; ./install $act 2>/dev/null
  cd $CASEDIR
  cd $SOURCE_ROOT_HYPRE; ./install $act 2>/dev/null
  cd $CASEDIR
}

function config_nek() {

# set defaults
: ${NEK_SOURCE_ROOT:="$HOME/Nek5000"}
: ${FC:="mpif90"}
: ${CC:="mpicc"}
: ${LD:=$FC}
: ${MPI:=1}
: ${PROFILING:=1}
: ${VISIT:=0}
: ${PPLIST:=""}
: ${NEK_DEPS_LOCAL_SRC:=0}

# first do some checks ...
if [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$PPLIST" == "?" ]; then
  echo "Options and arguments:"
  echo "  <casename>    name of usr file without extension"
  echo "  clean         cleanup build environment"
  echo "  distclean     cleanup build environment and remove tar files"
  echo "  -build-dep    build dependencies"
  echo ""
  echo "PPLIST config options:" 
  echo "  NOMPIIO       disable MPI-IO"
  echo "  XSMM          use libxsmm for mxm"
  echo "  CVODE         enable CVODE support"
  echo "  HYPRE         enable HYPRE support"
  echo "  VENDOR_BLAS   use VENDOR BLAS/LAPACK"
  echo "  DPROCMAP      use distributed processor mapping"
  echo "  PARRSB        use online RSB partitioner"
  return 11
fi


if [ "$VERMAKENEK" != "19.0.1" ]; then
  echo "FATAL ERROR: Your makenek is incompatible!"
  echo "Please upgrade to the version shipped the source"
  return 11
fi

rm -rf build.log
nargs=$#

CASEDIR=`pwd`
CASENAME=$1
NODEP=0
ONLYDEP=0
NEK_WARN=""

SOURCE_ROOT_BLAS="$NEK_SOURCE_ROOT/3rd_party/blasLapack"
export SOURCE_ROOT_BLAS
SOURCE_ROOT_GSLIB="$NEK_SOURCE_ROOT/3rd_party/gslib"
export SOURCE_ROOT_GSLIB
SOURCE_ROOT_PARRSB="$NEK_SOURCE_ROOT/3rd_party/parRSB"
export SOURCE_ROOT_PARRSB
SOURCE_ROOT_CVODE="$NEK_SOURCE_ROOT/3rd_party/cvode"
export SOURCE_ROOT_CVODE
SOURCE_ROOT_HYPRE="$NEK_SOURCE_ROOT/3rd_party/hypre"
export SOURCE_ROOT_HYPRE
export NEK_DEPS_LOCAL_SRC

NODEP=1
if [ "$1" == "-build-dep" ]; then
  ONLYDEP=1
  nargs=0
fi 

if [[ $1 = "clean" || $1 = "distclean" ]]; then
  while true; do
    read -p "do you want to cleanup all 3rd party dependencies too? [N]" yn
    case $yn in
      [Yy]* )
        clean_3rdparty "$1";
        if [ -f makefile ]; then make clean; fi
        break ;;
      * ) 
        if [ -f makefile ]; then make clean; fi
        break ;;
    esac
  done
  rm -rf build.log
  rm -rf makefile
  rm -rf .state
  return 0
fi

#if [ -f makefile ] && [ $ONLYDEP -ne 1 ]; then
#  while true; do
#    echo "makefile exists already!"
#    read -p "Are you sure you want to continue? [N]" yn
#    case $yn in
#      [Yy]* )
#        break ;;
#      * ) 
#        return 0
#        break ;;
#    esac
#  done
#fi

if [ ! -f $NEK_SOURCE_ROOT/core/makefile.template ]; then
  echo "FATAL ERROR: Cannot find $NEK_SOURCE_ROOT/core/makefile.template!"
  return 11
fi

# test compiler
if ! which `echo $FC | awk '{print $1}'` 1>/dev/null ; then
  echo "FATAL ERROR: Cannot find `echo $FC | awk '{print $1}'`!"
  return 11
fi
if ! which `echo $CC | awk '{print $1}'` 1>/dev/null ; then 
  echo "FATAL ERROR: Cannot find `echo $CC | awk '{print $1}'`!"
  return 11
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
  *pgf*)       FCPP="-Mpreprocess"
               FR8="-r8"
               FF77="-Mfixed"
               MODFLAG='-module $(OBJDIR)'
               ;;
  *flang*)
    if is_llvm_flang "$FCcomp"; then
        FCPP="-cpp"
        FR8="-fdefault-real-8 -fdefault-double-8"
        FF77="-w -ffixed-form"
        MODFLAG='-J$(OBJDIR)'
    else
        # CLASSIC FLANG (Legacy/PGI-based)
        FCPP="-Mpreprocess"
        FR8="-r8"
        FF77="-Mfixed"
        MODFLAG='-module $(OBJDIR)'
    fi
    ;;
  *nvfortran*) FCPP="-Mpreprocess"
               FR8="-r8"
               FF77="-Mfixed"
               MODFLAG='-module $(OBJDIR)'
               ;;
  *gfortran*)  FCPP="-cpp"
               FR8="-fdefault-real-8 -fdefault-double-8"
               FF77="-w -std=legacy -ffixed-form"
               MODFLAG='-J$(OBJDIR)'
               ;;
  *ftn*)       FCPP="-eZ"
               FR8="-sreal64"
               FF77="-fixed"
               MODFLAG='-J $(OBJDIR)'
               ;;
  *frt*)       FCPP="-Cpp"
               FR8="-CcdRR"
               FF77="-fixed"
               MODFLAG='-M $(OBJDIR)'
               ;;
  *armflang*)  FCPP="-cpp"
               FR8="-r8"
               FF77="-std=legacy -ffixed-form"
               MODFLAG='-J$(OBJDIR)'
               ;;
  *ifort*)     FCPP="-fpp"
               FR8="-r8 -fpconstant"
               FF77="-fixed"
               MODFLAG='-module $(OBJDIR)'
               ;;
  *ifx*)       FCPP="-fpp"
               FR8="-r8 -fpconstant"
               FF77="-fixed"
               MODFLAG='-module $(OBJDIR)'
               ;;
  *xlf*)       FCPP="-qsuffix=cpp=f"
               FR8="-qrealsize=8 -qdpc=e"
               FF77="-qfixed"
               PPPO="-WF,"
               MODFLAG='-qmoddir=$(OBJDIR)'
               ;;
  *)  echo "ERROR: Cannot find a supported compiler!"
      echo ""
      read;;
esac

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

# Check if ld option is available
cat > _test.f << _ACEOF
      program test
      end
_ACEOF

$FCcomp -Wl,--allow-multiple-definition -o _test _test.f >/dev/null 2>&1 || true
if [ -f _test ]; then
  LDFLAGS+=" -Wl,--allow-multiple-definition"
  \rm _test 2>/dev/null
fi

if [ "$(uname)" = "Darwin" ]; then
  $FCcomp -Wl,-ld_classic -o _test _test.f >/dev/null 2>&1 || true
  if [ -f _test ]; then
    LDFLAGS+=" -Wl,-ld_classic"
    \rm _test 2>/dev/null
  fi
fi

\rm _test.f 2>/dev/null

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
  PPLIST="${PPLIST} TIMER"
fi

PPLIST="${PPLIST}"

# add user defined preprocessor symbols
MPIIO=1
PARRSB=0
CVODE=0
HYPRE=0
MXM_USER="mxm_std.o"
VENDOR_BLAS=0 

for i in $PPLIST
do
   if [ "$i" == "NOMPIIO" ]; then
     MPIIO=0
   fi
   if [ "$i" == "PARRSB" ]; then
     PARRSB=1
     PARRSB_OPT=" DESTDIR=${SOURCE_ROOT_PARRSB} UNDERSCORE=$UNDERSCORE MPI=$MPI BLAS=0"
     PARRSB_OPT+=" GSLIBPATH=${SOURCE_ROOT_GSLIB}"
     USR_LFLAGS+=" -L${SOURCE_ROOT_PARRSB}/lib -lparRSB"
     PARRSB_INC="${SOURCE_ROOT_PARRSB}/include"
   fi
   if [ "$i" == "CVODE" ]; then 
     CVODE=1
     if [ $MPI -eq 1 ]; then
       USR_LFLAGS+=" -L$SOURCE_ROOT_CVODE/lib -lsundials_fcvode -lsundials_cvode"
       USR_LFLAGS+=" -lsundials_fnvecparallel -lsundials_nvecparallel" 
     else
       USR_LFLAGS+=" -L$SOURCE_ROOT_CVODE/lib -lsundials_fcvode -lsundials_cvode"
       USR_LFLAGS+=" -lsundials_fnvecserial -lsundials_nvecserial" 
     fi
   fi
   if [ "$i" == "HYPRE" ]; then
     HYPRE=1
     USR_LFLAGS+=" -L$SOURCE_ROOT_HYPRE/lib -lHYPRE"
     HYPRE_INC="$SOURCE_ROOT_HYPRE/include"
   fi
   if [ "$i" == "VENDOR_BLAS" ]; then
     VENDOR_BLAS=1
   fi
done

if [ $VENDOR_BLAS -ne 1 ]; then
  USR_LFLAGS+=" -L$SOURCE_ROOT_BLAS -lblasLapack"
  PARRSB_OPT+=" BLASDIR=${SOURCE_ROOT_BLAS}"
  PARRSB_OPT+=" BLASFLAGS=-lblasLapack"
fi

# gslib build options
GSLIB_OPT+=" DESTDIR=${SOURCE_ROOT_GSLIB}"
GSLIB_OPT+=" MPI=$MPI MPIIO=$MPIIO"
GSLIB_OPT+=" ADDUS=$UNDERSCORE USREXIT=0 BLAS=2"
GSLIB_IFLAGS=" -I$SOURCE_ROOT_GSLIB/include"
USR_LFLAGS+=" -L$SOURCE_ROOT_GSLIB/lib -lgs"

if ! which tar >/dev/null; then
  echo "FATAL ERROR: Cannot find tar!"
  return 11
fi
if ! which make >/dev/null; then
  echo "FATAL ERROR: Cannot find make!"
  return 11
fi
if [ $HYPRE -ne 0 ]; then
  if ! which cmake >/dev/null; then
    echo "FATAL ERROR: Cannot find cmake!"
    return 11
  fi
fi

export FC
export FCPP
export FR8
export FF77
export CC
export MPI
export MPIIO
export UNDERSCORE
export GSLIB_OPT
export PARRSB_OPT

if [ $ONLYDEP -eq 1 ]; then
  build_3rdparty
  return 0 
fi 

if [ $nargs -eq 0 ]; then
   NCASEFILES=$(find *.usr ! -name 'SIZE.usr' 2>/dev/null | wc -l)
   if [ $NCASEFILES == "0" ]; then
     true
   elif [ $NCASEFILES != "1" ]; then
     echo "FATAL ERROR: Found multiple .usr files!"
     echo "Please use makenek [.usr filename]"
     find *.usr ! -name 'SIZE.usr'
     return 11
   else
     CASENAME=$(find *.usr ! -name 'SIZE.usr' | sed 's/\.usr$//')
   fi
fi

if [ ! -f ${CASENAME}'.usr' ]; then
   echo "FATAL ERROR: Cannot find $1.usr!"
   return 11
fi

: ${INCLUDES:=""}
OPT_INCDIR="./ ${INCLUDES}"

printf "generating makefile ... "
rm -rf makefile 2>/dev/null

sed -e "s:^FC[ ]*=.*:FC=$FC:" \
-e "s:^CC[ ]*=.*:CC=$CC:" \
-e "s:^LD[ ]*=.*:LD=$LD:" \
-e "s:^OPT_INCDIR[ ]*=.*:OPT_INCDIR=$OPT_INCDIR:" \
-e "s:^FFLAGS[ ]*+=.*:FFLAGS+=$FCPP $FR8 $FF77 $FFLAGS:" \
-e "s:^MODFLAG[ ]*=.*:MODFLAG=$MODFLAG:" \
-e "s:^CFLAGS[ ]*+=.*:CFLAGS+=$CFLAGS:" \
-e "s:^LDFLAGS[ ]*+=.*:LDFLAGS+=$LDFLAGS:" \
-e "s/^PPPO[ ]*=.*/PPPO=$PPPO/" \
-e "s/^PPS[ ]*=.*/PPS=$PPLIST/" \
-e "s:^MXM[ ]*=.*:MXM=$MXM_USER:" \
-e "s:^BLAS[ ]*=.*:BLAS=$BLAS:" \
-e "s/^MPI[ ]*=.*/MPI=$MPI/" \
-e "s/^CVODE[ ]*=.*/CVODE=$CVODE/" \
-e "s:^GSLIB_IFLAGS[ ]*=.*:GSLIB_IFLAGS\:=${GSLIB_IFLAGS}:" \
-e "s/^VISIT[ ]*=.*/VISIT=${VISIT}/" \
-e "s:^USR[ ]*=.*:USR=$USR:" \
-e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS=$USR_LFLAGS:" \
-e "s:^S[ ]*=.*:S=${NEK_SOURCE_ROOT}:" \
$NEK_SOURCE_ROOT/core/makefile.template >.makefile

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

if [ $HYPRE -ne 0 ]; then
  sed -e "s:^HYPRE_IFLAGS[ ]*=.*:HYPRE_IFLAGS\:=-I${HYPRE_INC}:" \
      -e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS=$USR_LFLAGS:" \
      .makefile > .makefile.tmp

  mv .makefile.tmp .makefile
fi

if [ $PARRSB -ne 0 ]; then
  sed -e "s:^PARRSB_IFLAGS[ ]*=.*:PARRSB_IFLAGS\:=-I${PARRSB_INC}:" \
      -e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS=$USR_LFLAGS:" \
      .makefile > .makefile.tmp

  mv .makefile.tmp .makefile
fi


if [ "$USR" != "" ]; then
  echo "###########################################################" >> .makefile
  echo "include makefile_usr.inc" >> .makefile
fi

if [ -f .makefile ]; then
  sed -e "1i\\
### makefile automatically created by makenek `date +"%m/%d/%Y %T"` ###" \
-e "s:${NEK_SOURCE_ROOT}:\$S:g" \
-e "s:^S[ ]*=.*:S\:=${NEK_SOURCE_ROOT}:" \
-e "s:^CASEDIR[ ]*=.*:CASEDIR\:=${CASEDIR}:" \
-e "s:^CASENAME[ ]*=.*:CASENAME\:=${CASENAME}:" .makefile > makefile 
  printf "done\n"
else
  printf "failed!\n"
  return 1 
fi
\rm .makefile 2>/dev/null

# backward support SIZE and SIZEu
# Transition aid: legacy SIZE -> canonical SIZE.usr.
#   - lone SIZE            : auto-migrate (mv SIZE SIZE.usr).
#   - both, identical      : drop the duplicate SIZE (warn, no prompt).
#   - both, differ, TTY    : prompt y/N (y=overwrite SIZE.usr with SIZE; N=keep).
#   - both, differ, no TTY : FATAL, ask the user to resolve (CI/scripted builds
#     never silently pick one -- a leftover case-dir SIZE also SHADOWS core/SIZE
#     at compile time, so an unresolved SIZE must not remain).
if [ -f SIZE ]; then
  if [ -f SIZE.usr ]; then
    if cmp -s SIZE SIZE.usr; then
      rm -f SIZE
      NEK_WARN+="WARNING: removed duplicate legacy SIZE (identical to SIZE.usr)\n"
    else
      ans=N
      if [ -t 0 ]; then
        printf "\n"
        printf "SIZE and SIZE.usr both exist and DIFFER.\n"
        printf "  SIZE.usr is canonical; SIZE is legacy.\n"
        printf "  y = overwrite SIZE.usr with legacy SIZE\n"
        printf "  N = keep SIZE.usr (default)\n"
        read -r -p "Overwrite SIZE.usr with legacy SIZE? [y/N] " ans </dev/tty
      fi
      case "$ans" in
        y|Y)
          mv -f SIZE SIZE.usr
          NEK_WARN+="WARNING: SIZE.usr overwritten by legacy SIZE (user choice)\n"
          ;;
        *)
          echo ""
          echo "FATAL ERROR: SIZE and SIZE.usr both exist and differ."
          echo "  Resolve the conflict before continuing. Options:"
          echo "    - keep SIZE.usr : rm SIZE      (or mv SIZE SIZE.bak)"
          echo "    - use SIZE      : mv SIZE SIZE.usr"
          echo "    - Merge changes : port changes from SIZE into SIZE.usr"
          return 11
          ;;
      esac
    fi
  else
    NEK_WARN+="WARNING: SIZE is renamed to SIZE.usr\n"
    mv SIZE SIZE.usr
  fi
fi
if [ ! -f SIZEu.usr ] && [ -f SIZEu ]; then
  NEK_WARN+="WARNING: SIZEu is renamed to SIZE.usru\n"
  mv SIZEu SIZE.usru
fi

# tweak SIZE.usr file
if [ -f SIZE.usru ]; then
   if [ -f SIZE.usr ]; then
      echo "FATAL ERROR: found SIZE.usru and SIZE.usr!"
      return 11
   fi
   NEK_WARN+="WARNING: SIZE.usru will be deprecated in a future release!\n"
   mv -v SIZE.usru SIZE.usr
fi

if [ ! -f SIZE.usr ]; then
   echo "FATAL ERROR: Cannot find SIZE.usr!"
   return 11
fi

if ! cat SIZE.usr | grep -qi 'lelr' ; then
   echo >>SIZE.usr
   echo 'c automatically added by makenek' >>SIZE.usr
   echo '      integer lelr' >>SIZE.usr
   echo '      parameter (lelr=lelt) ! max number of local elements per restart file' >>SIZE.usr
fi
if ! cat SIZE.usr | grep -qi 'ldimt_proj' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer ldimt_proj' >>SIZE.usr
     echo '      parameter(ldimt_proj=1) ! max auxiliary fields residual projection ' >>SIZE.usr
fi

# tweak leagacy SIZE.usr file
if ! cat SIZE.usr | grep -q 'SIZE.inc' ; then
  NEK_WARN+="WARNING: Not include SIZE.inc in SIZE.usr is deprecated and may be unsupported later.\n"
  if ! cat SIZE.usr | grep -qi 'optlevel' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer optlevel,loglevel' >>SIZE.usr
     echo '      common /lolevels/ optlevel,loglevel' >>SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'lxo' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer lxo' >> SIZE.usr
     echo '      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)' >>SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'ax1' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer ax1,ay1,az1,ax2,ay2,az2' >> SIZE.usr
     echo '      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! running averages' >> SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'lys=lxs' ; then
     cat SIZE.usr | grep -iv lxs > SIZE.usr.x; mv SIZE.usr.x SIZE.usr  # Clean existing SIZE.usr file of old version
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer lxs,lys,lzs' >> SIZE.usr
     echo '      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Preconditioner' >> SIZE.usr

  fi
  if ! cat SIZE.usr | grep -qi 'lcvx1' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer lcvx1,lcvy1,lcvz1,lcvelt' >> SIZE.usr
     echo '      parameter (lcvx1=1,lcvy1=1,lcvz1=1,lcvelt=1) ! cvode arrays' >> SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'lfdm' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer lfdm' >> SIZE.usr
     echo '      parameter (lfdm=0)  ! == 1 for fast diagonalization method' >> SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'nsessmax' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer nsessmax' >> SIZE.usr
     echo '      parameter (nsessmax=1)  ! max sessions to NEKNEK' >> SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'nmaxl_nn' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer nmaxl_nn' >> SIZE.usr
     echo '      parameter (nmaxl_nn=' >> SIZE.usr
     echo '     $          min(1+(nsessmax-1)*2*ldim*lxz*lelt,2*ldim*lxz*lelt))' >>SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'nfldmax_nn' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer nfldmax_nn' >> SIZE.usr
     echo '      parameter (nfldmax_nn=' >> SIZE.usr
     echo '     $          min(1+(nsessmax-1)*(ldim+1+ldimt),ldim+1+ldimt))' >>SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'nio' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer nio' >> SIZE.usr
     echo '      common/IOFLAG/nio  ! for logfile verbosity control' >> SIZE.usr
  fi
  if ! cat SIZE.usr | grep -qi 'lhref' ; then
     echo >>SIZE.usr
     echo 'c automatically added by makenek' >>SIZE.usr
     echo '      integer lhref' >> SIZE.usr
     echo '      parameter (lhref=10)' >> SIZE.usr
  fi
fi

declare -a varlist=(
"NEK_SOURCE_ROOT"
"FC"
"FFLAGS"
"FCPP"
"FR8"
"FF77"
"CC"
"CFLAGS"
"PPLIST"
)
if [ -f .state ]; then
  for i in "${varlist[@]}"
  do
    if ! cat .state | grep -q "$i=${!i}" ; then
      echo "configuration change detected, preparing fresh build"
      rm -rf .state
      while true; do
        read -p "do you want to rebuild all 3rd party dependencies too? [N]" yn
        case $yn in
          [Yy]* )
            clean_3rdparty "clean";
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

if [ $NODEP -ne 1 ]; then
  build_3rdparty
fi

# Collet NEK_WARN to obj/.nekwarn for later echo
mkdir -p "$CASEDIR/obj"
\rm -f "$CASEDIR/obj/.nekwarn"
if [ -n "$NEK_WARN" ]; then
   printf '%b' "$NEK_WARN" > "$CASEDIR/obj/.nekwarn"
   # Also echo now (before make runs) so warnings are visible even if the
   # compilation fails; on success they are shown again before the banner.
   echo ""
   printf '%b' "$NEK_WARN"
   echo ""
fi

return 0
}
