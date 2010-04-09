# This include file is used by the makenek script
# to create a makefile using a makefile template 
# ------------------------------------------------

# first do some checks ...
if [ $# -eq 0 ]; then
  echo ""
  echo "makenek - automatic build tool for nek5000"
  echo "usage: makenek [.usr filename | clean]"
  echo ""
  exit 1
fi

CASENAME=$1
CASEDIR=`pwd`
APATH_SRC=`cd $SOURCE_ROOT; pwd`
SOURCE_ROOT=$APATH_SRC

if [ $1 = "clean" ]; then
  make clean
  exit 0
fi

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

which `echo $CC | awk '{print $1}'` 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $CC!"
  exit 1
fi

\rm test_cc.o 2>/dev/null
# basic compiler test
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
rm -f $SOURCE_ROOT/SIZE      2>/dev/null
rm -f $SOURCE_ROOT/subuser.f 2>/dev/null
rm -f ./obj/subuser.o        2>/dev/null
\mv ./obj/*.o .              2>/dev/null

# Check if the compiler adds an underscore to external functions
IFUNDERSCORE=false
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call byte_write
      end
_ACEOF
$F77 -c test_underscore.f 2>&1 >/dev/null 
nm test_underscore.o | grep byte_write_ 1>/dev/null
if [ $? -eq 0 ] 
then
  IFUNDERSCORE=true
fi
\rm test_underscore.* 2>/dev/null

# assign compiler flag to promote REAL to 8 bytes 
F77comp=$F77
echo $F77comp | grep -i mpi >/dev/null
if [ $? -eq 0 ]; then   
 $F77 -showme  &>/dev/null
 if [ $? -eq 0 ]; then
    F77comp=`$F77 -showme | awk '{print $1}' | awk -F/ '{print $NF}'`
 else
    F77comp=`$F77 -show | awk '{print $1}' | awk -F/ '{print $NF}'`
 fi
else
 F77comp=`echo $F77 | awk '{print $1}'`
fi

case $F77comp in
  pgf*)        P="-r8 -Mpreprocess"
               ;;
  gfortran)    P="-fdefault-real-8 -x f77-cpp-input"
               ;;
  ifort)       P="-r8 -fpconstant -fpp"
               ;;
  pathf*)      P="-r8 -cpp -fno-second-underscore"
               ;;
  *xlf*)       P="-qrealsize=8 -qdpc=e -qsuffix=cpp=f"
               PPPO="-WF,"
               ;;
  ftn)         P="-r8 -Mpreprocess"
               ;;
  sun*)        P="-r8const -xtypemap=real:64 -fpp"
               ;;
  open*)       P="-r8 -cpp -fno-second-underscore"
               ;;
  *)  echo "WARNING: Cannot detect compiler specfic flags
      echo "          - to promote REAL to 8 bytes"
      echo "          - to invoke preprocessor first
      echo "Please edit the makefile and specify the requested compiler flags in the P variable!"
      echo ""
      P=""
      echo "Press return to finish the makefile generation ..."
      read;;
esac
export PPPO

# set logical flags
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

IFMOAB=false

# add user defined preprocessor symbols
for i in $PPLIST
do
   P="$P $PPPO-D$i"
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
     fi
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
    cat SIZE | grep -i 'cv_maxl' >/dev/null
    if [ $? -ne 0 ]; then
      cp SIZE .SIZE
      echo >>.SIZE
      echo 'c CVODE settings automatically added by makenek' >>.SIZE
      echo '      integer cv_maxl' >> .SIZE
      echo '      real cv_delt,cv_sigs' >> .SIZE
      echo '      parameter(cv_maxl = 20  ) ! Size of Krylov Space'  >>.SIZE
      echo '      parameter(cv_delt = 0.1 ) ! linear convergence factor '  >>.SIZE
      echo '      parameter(cv_sigs = 0.01) ! Jacobian DQ perturbation scaling factor'  >>.SIZE
      mv .SIZE SIZE
    fi
   fi
done


# set preprocessor symbol to print AMG files
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

# set preprocessor symbol for timing
if [ "$IFTIMING" == "false" -o "$IFTIMING" == "no" ]; then
    P="$P $PPPO-DNOTIMER"
fi

if [ "$IFAMG" == "true" -o "$IFAMG" == "yes" ]; then
  IFAMG=true
else
  # default
  IFAMG=false
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

# Generate makefile from using the makefile template and 
# the specified parameters in makenek

echo "generating makefile ..."
rm -rf makefile 2>/dev/null

# MXM-handling
if [ "$MXM_USER" != "" ]; then
    echo "ABORT: makenek version is too old! Update with latest from source!"
    exit 1 
fi
MXM_USER="mxm_std.o blas.o"
echo $PPLIST | grep 'BG_MXM' >/dev/null 
if [ $? -eq 0 ]; then
   MXM_USER="mxm_std.o bg_aligned3.o bg_mxm44.o bg_mxm44_uneven.o bg_mxm3.o blas.o" 
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


sed -e "s:^F77[ ]*=.*:F77\:=$F77:" \
-e "s:^CC[ ]*=.*:CC\:=$CC:" \
-e "s:^G[ ]*=.*:G\:=$G:" \
-e "s:^OPT_FLAGS[ ]*=.*:OPT_FLAGS\:=$OPT_FLAGS:" \
-e "s/^P[ ]*=.*/P:=$P/" \
-e "s/^L0[ ]*=.*/L0=$L0/" \
-e "s/^L2[ ]*=.*/L2=$L2/" \
-e "s/^L3[ ]*=.*/L3=$L3/" \
-e "s/^L4[ ]*=.*/L4=$L4/" \
-e "s:^MXM[ ]*=.*:MXM=$MXM_USER:" \
-e "s/^IFAMG[ ]*=.*/IFAMG:=$IFAMG/" \
-e "s/^IFAMG_DUMP[ ]*=.*/IFAMG_DUMP:=$IFAMG_DUMP/" \
-e "s/^IFMPI[ ]*=.*/IFMPI:=$IFMPI/" \
-e "s/^IFMOAB[ ]*=.*/IFMOAB:=$IFMOAB/" \
-e "s:^MOAB_DIR[ ]*=.*:MOAB_DIR\:=${MOAB_DIR}:" \
-e "s/^IFUNDERSCORE[ ]*=.*/IFUNDERSCORE:=$IFUNDERSCORE/" \
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

#mv $SOURCE_ROOT/HEADER $SOURCE_ROOT/.HEADER 
#which svnversion >/dev/null 2>&1 
#status=`echo $?`
#if [ "$status" == "0" ]; then
#  rev_number=$( svnversion ${SOURCE_ROOT} )
#  rev_number=$( printf %-11s ${rev_number} )
#  sed "s/SVN r.........../SVN r$rev_number/" \
#    $SOURCE_ROOT/.HEADER > $SOURCE_ROOT/HEADER
#else
#  sed "s/SVN r.........../SVN r\?\?\?\?/" \
#    $SOURCE_ROOT/.HEADER > $SOURCE_ROOT/HEADER
#fi
#rm -rf $SOURCE_ROOT/.HEADER >/dev/null 2>&1

# add plugins
for i in $PLUGIN_LIST
do
  if [ ! -f makenek_plugin_$i ]; then
     echo "ABORT: cannot find plug-in config file" makenek_plugin_$i
     exit 1
  fi
  source makenek_plugin_$i

  plugindir="SOURCE_ROOT_$i"
  plugindir=`echo ${!plugindir}`

  if [ ! -f $plugindir/makenek.inc ]; then
     echo "ABORT: Cannot find plug-in include file file" $plugindir/makenek.inc
     exit 1
  fi
  source $plugindir/makenek.inc
  if [ "$ierr" != 0 ]; then
     echo "ABORT: plug-in config failed!"
     exit 1
  fi
  echo "found plug-in" $i
done

if [ -f SIZE ]; then 
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
else
  echo "FATAL ERROR: Cannot find SIZE"
  exit 1
fi
