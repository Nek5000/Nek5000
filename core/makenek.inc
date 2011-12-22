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
        call byte_write
      end
_ACEOF
$F77 -c test_underscore.f 2>&1 >/dev/null 
nm test_underscore.o | grep byte_write_ 1>/dev/null
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
  *gfortran*)   P="-fcray-pointer -fdefault-real-8 -x f77-cpp-input"
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
-e "s/^IFAMG[ ]*=.*/IFAMG:=$IFAMG/" \
-e "s/^IFAMG_DUMP[ ]*=.*/IFAMG_DUMP:=$IFAMG_DUMP/" \
-e "s/^IFMPI[ ]*=.*/IFMPI:=$IFMPI/" \
-e "s/^IFMOAB[ ]*=.*/IFMOAB:=$IFMOAB/" \
-e "s:^MOAB_DIR[ ]*=.*:MOAB_DIR\:=${MOAB_DIR}:" \
-e "s/^IFNEKNEK[ ]*=.*/IFNEKNEK:=$IFNEKNEK/" \
-e "s:^USR[ ]*=.*:USR\:=$USR:" \
-e "s:^USR_LFLAGS[ ]*=.*:USR_LFLAGS\:=$USR_LFLAGS:" \
-e "s:^S[ ]*=.*:S\:=${SOURCE_ROOT}:" $SOURCE_ROOT/makefile.template >.makefile

echo $G | grep '\-g' 1>/dev/null
if [ $? -eq 0 ]; then
  sed 's/-O[1-4]/-O0/g' .makefile > .makefile.tmp
  mv .makefile.tmp .makefile
  echo "Activate DEBUG mode"
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
     echo '      parameter(lpart = 1  ) ! max number of particles' >>SIZE
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

if [ $NOCOMPILE -eq 1 ]; then
  exit 0
fi 
