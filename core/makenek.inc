# This is include file is needed for the makenek script
# -----------------------------------------------------

# do some checks

if [ $# -eq 0 ]
then
  echo "missing argument!"
  echo "usage: makenek [.usr filename] [clean]"
  exit 1
fi

if [ ! -f $SOURCE_ROOT/makefile ]
then
  echo "$SOURCE_ROOT/makefile does not exist!"
  exit 1
fi

which $F77 1>/dev/null
if [ $? -ne 0 ]
then
  echo "Cannot find $F77! Check settings in makenek."
  exit 1
fi

which $CC 1>/dev/null
if [ $? -ne 0 ]
then
  echo "Cannot find $CC! Check settings in makenek."
  exit 1
fi

if [ $# -eq 1 ]
then
  if [ $1 = "clean" ]
  then
    make -f $SOURCE_ROOT/makefile clean
    exit 0
  fi
  
  if [ ! -f ${1}'.usr' ]
  then
     echo "Cannot find $1.usr!"
     exit 1
  fi
fi

rm -f nek5000 2>/dev/null
rm -f $SOURCE_ROOT/SIZE 2>/dev/null
rm -f $SOURCE_ROOT/subuser.f 2>/dev/null

\mv ./obj/*.o . 2>/dev/null
ln -sf `pwd`/SIZEu    $SOURCE_ROOT/SIZE
ln -sf `pwd`/$1.usr   $SOURCE_ROOT/subuser.f

# check if the compiler adds an underscore to external functions
cat > tmp.f << _ACEOF
      subroutine tmp
        call byte_write
      end
_ACEOF

$F77 -c tmp.f 1>/dev/null 
nm tmp.o | grep byte_write_ 1>/dev/null
if [ $? -eq 0 ] 
then
  UNDERSCORE=" -DUNDERSCORE"
fi
rm -f tmp.o 2>/dev/null

# check flag for casting to 8 byte reals
cat > tmp.f << _ACEOF
      program tmp
        if ((1.+1.e-12).ne.1.) write(*,'(A)') 'wdsize=8'
      end
_ACEOF

$F77 -r8 tmp.f 2>$1>/dev/null 

if [ -f "a.out" ]; then status=`./a.out`; fi
 
if [ "$status" == "wdsize=8" ]
then
  P=" -r8"
else
  xlf -qrealsize=8 tmp.f 2>$1>/dev/null
  if [ -f "a.out" ]; then status=`./a.out`; fi
  if [ "$status" == "wdsize=8" ] 
  then
    P=" -qrealsize=8"
  else
   echo "cannot detect compiler flag for 8 byte reals"
   exit 1
  fi
fi
rm tmp.f a.out 2>/dev/null

export F77
export CC
export G
export P
export IFMPI
export UNDERSCORE
export SOURCE_ROOT
export USR
export USR_LIB
