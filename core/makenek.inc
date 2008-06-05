# This include file is used by the makenek script
# to create a makefile using a makefile template 
# ------------------------------------------------

# first do some checks ...
if [ $# -eq 0 ]; then
  echo "ERROR: missing argument!"
  echo "usage: makenek [.usr filename] [clean]"
  exit 1
fi

if [ $# -eq 1 ]; then
  if [ $1 = "clean" ]; then
    make -f $SOURCE_ROOT/makefile.template clean
    exit 0
  fi
  if [ ! -f ${1}'.usr' ]; then
     echo "FATAL ERROR: Cannot find $1.usr!"
     exit 1
  fi
fi

if [ ! -f $SOURCE_ROOT/makefile.template ]; then
  echo "FATAL ERROR: Cannot find $SOURCE_ROOT/makefile.template!"
  exit 1
fi

which $F77 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $F77!"
  exit 1
fi

# basic compiler test
cat > tmp.f << _ACEOF
      program test
      end
_ACEOF
rm a.out 2>/dev/null
$F77 tmp.f 2>&1 > /dev/null
if [ ! -f a.out ]; then
  echo "FATAL ERROR: Basic compiler test for $F77 failed!"
  exit 1
fi
rm tmp.f a.out 2>/dev/null

which $CC 1>/dev/null
if [ $? -ne 0 ]; then
  echo "FATAL ERROR: Cannot find $CC!"
  exit 1
fi

# basic compiler test
cat > tmp.c << _ACEOF
      int main(){}
_ACEOF
rm a.out 2>/dev/null
$CC tmp.c
if [ ! -f a.out ]; then
  echo "FATAL ERROR: Basic compiler test for $CC failed!"
  exit 1
fi
rm tmp.c a.out 2>/dev/null


rm -f nek5000 2>/dev/null
rm -f $SOURCE_ROOT/SIZE 2>/dev/null
rm -f $SOURCE_ROOT/subuser.f 2>/dev/null

\mv ./obj/*.o . 2>/dev/null
cp  `pwd`/SIZEu    $SOURCE_ROOT/SIZE
cp  `pwd`/$1.usr   $SOURCE_ROOT/subuser.f


# Check if the compiler adds an underscore to external functions
cat > tmp.f << _ACEOF
      subroutine underscore_test
        call byte_write
      end
_ACEOF

IFUNDERSCORE=false
$F77 -c tmp.f 1>/dev/null 
nm tmp.o | grep byte_write_ 1>/dev/null
if [ $? -eq 0 ] 
then
  IFUNDERSCORE=true
fi
rm -f tmp.o 2>/dev/null


# Check compiler flag for casting to 8 byte reals
cat > tmp.f << _ACEOF
      program tmp
        if ((1.+1.e-12).ne.1.) write(*,'(A)') 'wdsize=8'
      end
_ACEOF

$F77 -r8 tmp.f 2>&1 >/dev/null 
if [ -f "a.out" ]; then status=`./a.out`; fi
if [ "$status" == "wdsize=8" ]
then
  P="-r8"
else
  xlf -qrealsize=8 tmp.f 2>&1 >/dev/null
  if [ -f "a.out" ]; then status=`./a.out`; fi
  if [ "$status" == "wdsize=8" ] 
  then
    P="-qrealsize=8"
  else
   echo "WARNING: Cannot detect compiler flag to promote REAL to 8 bytes!"
   echo "         Please specify the compiler option using the P flag in the makefile"
   echo "Press return to continue"
   read
  fi
fi
rm tmp.f a.out 2>/dev/null

# set logical flags
if [ "$IFMPI" == "false" -o "$IFMPI" == "no" ]; then
  IFMPI=false
else
  # default
  IFMPI=true
fi

if [ "$IFAMG" == "true" -o "$IFAMG" == "yes" ]; then
  IFAMG=true
else
  # default
  IFAMG=false
fi

# Generate makefile from using the makefile template and 
# the specified parameters in makenek

echo "try to generate makefile ..."

sed -e "s:^F77[ ]*=.*:F77\:=$F77:" \
-e "s:^CC[ ]*=.*:CC\:=$CC:" \
-e "s/^G[ ]*=.*/G:=$G/1" \
-e "s/^P[ ]*=.*/P:=$P/" \
-e "s/^IFAMG[ ]*=.*/IFAMG:=$IFAMG/" \
-e "s/^IFMPI[ ]*=.*/IFMPI:=$IFMPI/" \
-e "s/^IFUNDERSCORE[ ]*=.*/IFUNDERSCORE:=$IFUNDERSCORE/" \
-e "s:^USR[ ]*=.*:USR\:=$USR:" \
-e "s:^USR_LIB[ ]*=.*:USR_LIB\:=$USR_LIB:" \
-e "s:^S[ ]*=.*:S\:=${SOURCE_ROOT}:" $SOURCE_ROOT/makefile.template >.makefile

echo $G | grep '\-g' 1>/dev/null
if [ $? -eq 0 ]; then
  sed 's/-O[1-4]/-O0/g' .makefile > .makefile.tmp
  mv .makefile.tmp .makefile
  echo "DEBUG mode activated - press return to continue"
  read
fi

if [ "$USR" != "" ]; then
  echo "###########################################################" >> makefile
  echo "include makefile_usr.inc" >> .makefile
fi

if [ -f .makefile ]; then
  sed "1i\### makefile automatically created by makenek `date +"%m/%d/%Y %T"` ###\n" \
      .makefile > makefile
  echo "makefile created"
else
  echo "FATAL ERROR: Makefile could not be created!" 
fi
\rm .makefile 2>/dev/null
