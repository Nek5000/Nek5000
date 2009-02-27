# This include file is used by the makenek script
# to create a makefile using a makefile template 
# ------------------------------------------------

# first do some checks ...
if [ $# -ne 1 ]; then
  echo "usage: makenek [.usr filename|clean]"
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
else
 F77comp=`echo $F77 | awk '{print $1}'`
fi

case $F77comp in
  pgf*)        P="-r8 -Mpreprocess"
               ;;
  gfortran)    P="-fdefault-real-8 -x f77-cpp-input"
               ;;
  ifort)       P="-r8 -fpp"
               ;;
  pathf*)      P="-r8 -cpp -fno-second-underscore"
               ;;
  *xlf*)       P="-qrealsize=8 -qsuffix=cpp=f"
               PPPO="-WF,"
               ;;
  ftn)         P="-r8 -Mpreprocess"
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



# set preprocssor symbol for MOAB support
if [ "$IFMOAB" == "true" -o "IFMOAB" == "yes" ]; then
   P="$P $PPPO-DMOAB"
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
else
   IFMOAB=false
fi

# set preprocssor symbol for timing
if [ "$IFTIMING" == "false" -o "$IFTIMING" == "no" ]; then
    P="$P $PPPO-DNOTIMER"
fi

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
    USR_LFLAGS="${USR_LFLAGS} $USR_LIB"
fi

# Generate makefile from using the makefile template and 
# the specified parameters in makenek

echo "try to generate makefile ..."
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
-e "s/^IFAMG[ ]*=.*/IFAMG:=$IFAMG/" \
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
  echo "DEBUG mode activated - press return to continue"
fi

if [ "$USR" != "" ]; then
  echo "###########################################################" >> makefile
  echo "include makefile_usr.inc" >> .makefile
fi


if [ "$MXM_USER" != "" ]; then
  mv .makefile .makefile.tmp
  sed -e "s:^MXM.*:MXM=$MXM_USER:" .makefile.tmp > .makefile 
fi

if [ -f .makefile ]; then
  sed "1i\### makefile automatically created by makenek `date +"%m/%d/%Y %T"` ###\n" \
      .makefile > makefile
  echo "makefile created"
else
  echo "FATAL ERROR: Makefile could not be created!" 
fi
\rm .makefile 2>/dev/null

mv $SOURCE_ROOT/HEADER $SOURCE_ROOT/.HEADER 
which svnversion >/dev/null 2>&1 
status=`echo $?`
if [ "$status" == "0" ]; then
  rev_number=$( svnversion ${SOURCE_ROOT} )
  rev_number=$( printf %-11s ${rev_number} )
  sed "s/SVN r.........../SVN r$rev_number/" \
    $SOURCE_ROOT/.HEADER > $SOURCE_ROOT/HEADER
else
  sed "s/SVN r.........../SVN r\?\?\?\?/" \
    $SOURCE_ROOT/.HEADER > $SOURCE_ROOT/HEADER
fi
rm -rf $SOURCE_ROOT/.HEADER >/dev/null 2>&1
