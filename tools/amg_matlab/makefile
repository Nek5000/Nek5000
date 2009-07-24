SHELL = /bin/sh

# USER SETTINGS ##################

# compilers
CC='pgcc'
F77='pgf77'
F90='pgf90'

# binary path
bin_nek_tools="/home/$(USER)/bin"

# tools to compile
MODULES = genMap genBox reatore2 postnek n2to3 prenek nekmerge








##################################
















all:
	 @export CC=${CC} F77=${F77} F90=${F90} bin_nek_tools=${bin_nek_tools};\
	 cp -v ./misc/* ${bin_nek_tools};   \
	 for i in ${MODULES} ; do           \
      if test -d $${i} ; then          \
        echo "----------------------"; \
        echo "Make $${i}..." ;         \
        echo "----------------------"; \
        cd $${i} ;                     \
        ${MAKE} ;                      \
        if test  $$? -ne 0 ; then      \
         exit 1;                       \
		  fi ; 									\
		  cd .. ;                        \
        echo "" ;                      \
      fi ;                             \
    done ;									   \
	echo "##############################################################";\
	echo "  TOOLS were installed in ${bin_nek_tools}                    ";\
	echo "##############################################################";

clean:
	@for i in ${MODULES} ; do         \
      if test -d $${i} ;  then       \
        cd $${i} ;                   \
        ${MAKE} clean ;              \
        cd .. ;                      \
      fi ;                           \
    done ;                           \
	 rm -f ${bin_nek_tools}/nekk      \
			 ${bin_nek_tools}/neks      \
 			 ${bin_nek_tools}/prex      \
 			 ${bin_nek_tools}/nekmerge  \
 			 ${bin_nek_tools}/reatore2  \
 			 ${bin_nek_tools}/genbox    \
  			 ${bin_nek_tools}/genmap    \
			 ${bin_nek_tools}/n2to3;    \
