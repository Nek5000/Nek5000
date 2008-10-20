# This include file is used by the makenek script
# to create a makefile using a makefile template 
# ------------------------------------------------

# first do some checks ...
if [ $# -ne 1 ]; then
  echo "usage: makenek <.usr filename | clean>"
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

which $CC 1>/dev/null
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

rm -f nek5000 2>/dev/null
rm -f $SOURCE_ROOT/SIZE 2>/dev/null
rm -f $SOURCE_ROOT/subuser.f 2>/dev/null

\mv ./obj/*.o . 2>/dev/null
if [ -f SIZEu ]; then 
  cp SIZEu    $SOURCE_ROOT/SIZE
else
  echo "FATAL ERROR: Cannot find SIZEu"
  exit 1
fi
cp $1.usr   $SOURCE_ROOT/subuser.f


# Check if the compiler adds an underscore to external functions
cat > test_underscore.f << _ACEOF
      subroutine underscore_test
        call byte_write
      end
_ACEOF

IFUNDERSCORE=false
$F77 -c test_underscore.f 2>&1 >/dev/null 
nm test_underscore.o | grep byte_write_ 1>/dev/null
if [ $? -eq 0 ] 
then
  IFUNDERSCORE=true
fi
\rm test_underscore.* 2>/dev/null

# assign compiler flag to promote REAL to 8 bytes 
F77comp=$F77
echo $F77comp | grep -i mpi
if [ $? -eq 0 ]; then
 F77comp=`$F77 -show | awk '{print $1}' | awk -F/ '{print $NF}'`
fi

case $F77comp in
  pgf77)       P="-r8";;
  pgf90)       P="-r8";;
  gfortran)    P="-fdefault-real-8";;
  ifort)       P="-r8";;
  bgxlf)       P="-qrealsize=8";;
  bgxlf_r)     P="-qrealsize=8";;
  blrts_xlf)   P="-qrealsize=8";;
  qk-pgf77)    P="-r8";;
  *)  echo "WARNING: Cannot detect compiler flag to promote REAL to 8 bytes!"
      echo "         Please specify the right compiler option using the P flag in the makefile"
      P=""
      echo "Press return to continue"
      read;;
esac

# set logical flags
if [ "$IFMPI" == "false" -o "$IFMPI" == "no" ]; then
  IFMPI=false
else
  # default
  IFMPI=true
fi
export IFMPI

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
-e "s:^G[ ]*=.*:G\:=$G:" \
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
  if [ "$P" = "" ]; then
    exit 1
  fi
else
  echo "FATAL ERROR: Makefile could not be created!" 
fi
\rm .makefile 2>/dev/null
