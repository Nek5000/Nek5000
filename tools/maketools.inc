#set -x
trap onexit EXIT

function onexit() {
  retval=$?
  if [ $retval -gt 0 ] && [ $retval -lt 11 ]; then
     printf "\nFATAL ERROR occured during compilation!\nCheck build.log for details\n"
  fi
}

: ${FC:="gfortran"}
: ${CC:="gcc"}
: ${bin_nek_tools:=`cd ../bin && pwd`}
: ${MAXNEL:=150000}

if [ $# -eq 0 ] || [ $1 == '-h' ] || [ $1 == '-help' ]; then
   echo "Usage: maketools [clean|all|tool(s)]"
   exit 0
fi

if [ "$FC" == "" ]; then
  echo "FATAL ERROR: Specify your Fortran compiler in maketools!"
  exit 11
fi

if [ "$CC" == "" ]; then
  echo "FATAL ERROR: Specify your C compiler in maketools!"
  exit 11
fi

which `echo $FC | awk '{print $1}'` 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $FC!"
  exit 11
fi

which `echo $CC | awk '{print $1}'` 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $CC!"
  exit 11
fi

if [ ! -d $bin_nek_tools  ]; then
  echo "FATAL ERROR: install path $bin_nek_tools does not exist!"
  exit 11
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

PPPO=""
case $FCcomp in
  *pgf*)       R8="-r8"
               CPPF="-Mpreprocess"
               BIGMEM="-mcmodel=medium"
               ;;
  *gfortran*)  R8="-fdefault-real-8 -fdefault-double-8"
               CPPF="-std=legacy -cpp"
               BIGMEM="-mcmodel=medium"
               ;;
  *ifort*)     R8="-r8"
               CPPF="-fpp"
               BIGMEM="-mcmodel=medium -shared-intel"
               ;;
  *xlf*)       R8="-qrealsize=8"
               CPPF="-qsuffix=cpp=f"
               BIGMEM="-q64"
               PPPO="-WF,"
               ;;
  *)           echo "FATAL ERROR: Cannot find a supported compiler!"
               exit 11
               ;;
esac

# Check if the compiler adds an underscore to external functions
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call byte_write
      end
_ACEOF

$FC -c test_underscore.f 2>&1 >/dev/null 
nm test_underscore.o | grep byte_write_ 1>/dev/null
if [ $? -eq 0 ]; then 
  US="-DUNDERSCORE"
fi
\rm test_underscore.* 2>/dev/null

# Test BIGMEM support
cat > _test.f << _ACEOF
      program test
      end
_ACEOF
$FC $BIGMEM -o _test _test.f >/dev/null 2>&1 || true
\rm _test.f 2>/dev/null
if [ -f _test ]; then
  \rm _test 2>/dev/null
else
 BIGMEM=""
fi

FFLAGS+=" -O2 -DMAXNEL=$MAXNEL"
FFLAGS=`echo $FFLAGS | sed -e "s/-D/$PPPO-D/g"`

export FC=`which $FC`
export NEK_FFLAGS="$FFLAGS"
export FFLAGS_
export CC=`which $CC`
export NEK_CFLAGS="$CFLAGS"
export bin_nek_tools
export US
export R8
export CPPF
export BIGMEM
export LDFLAGS

#LIST=`ls -Cd */ | sed 's:\/::g'`
LIST="genmap gencon genbox n2to3 reatore2 nekmerge prenek postnek nekamg_setup gmsh2nek exo2nek cgns2nek"
export LIST

rm -rf build.log
if [ "$1" == "clean" ]; then
   for i in $LIST ; do
     make -C ./$i clean >/dev/null 2>&1
     rm -rf ./$i/build.log
     rm -rf ../bin/$i
   done 
else
   export DOALL=''
   if [ "$1" == "all" ]; then
      export DOALL=1
   fi
   if [ "$1" != "all" ]; then
      LIST=`echo $* | sed 's:\/::g'`
   fi
   set -o pipefail
   for TOOL in $LIST ; do
     export FFLAGS="${NEK_FFLAGS}" 
     export FFLAGS+=" ${CPPF}"
     export FFLAGS+=" ${BIGMEM}"
     export CFLAGS="${NEK_CFLAGS}"
     export CFLAGS+=" ${US}"
     export CFLAGS+=" ${BIGMEM}"
     if [ "$TOOL" == "exo2nek" ]; then
        export FFLAGS+=" ${R8}"
     fi
     if [ "$TOOL" == "gmsh2nek" ]; then
        export FFLAGS+=" ${R8}"
     fi	 
     if [ "$TOOL" == "n2to3" ]; then 
        export FFLAGS+=" ${R8}"
     fi
     if [ "$TOOL" == "genmap" ]; then
        export FFLAGS+=" ${R8}"
     fi
     if [ "$TOOL" == "gencon" ]; then
        export FFLAGS+=" ${R8}"
     fi
     if [ "$TOOL" == "genbox" ]; then
       export FFLAGS+=" ${R8}"
     fi
     if [ "$TOOL" == "prenek" ]; then
        export CFLAGS+=" -Dr8"
     fi
     printf "building $TOOL ... "
     if make -C ./$TOOL 2>&1 | tee ./$TOOL/build.log >/dev/null 2>&1; then
       printf "done\n"
     else
       printf "failed\n"
     fi
   done
   set +o pipefail 
fi
