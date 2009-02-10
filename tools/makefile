SHELL = /bin/sh

all:
	 @for i in $(MODULES) ; do               \
      if test -d $${i} ; then               \
        echo "----------------------";      \
        echo "Make $${i}..." ;              \
        echo "----------------------";      \
        cd $${i} ;                          \
        export FFLAGS_IN=  ;                \
        export CFLAGS_IN=${US} ;            \
        if [ "$$i" == "genMap" ]; then      \
           export FFLAGS_IN=${R8} ;         \
        fi ;                                \
        if [ "$$i" == "prenek" ]; then      \
           export FFLAGS_IN=${R8} ;         \
           export CFLAGS_IN="$(US) -Dr8";   \
        fi ;                                \
        ${MAKE} ;                           \
        if test  $$? -ne 0 ; then           \
          exit 1 ;                          \
        fi ;                                \
		  cd .. ;                             \
        echo "" ;                           \
      fi ;                                  \
    done ;									        \
	echo "##############################################################";\
	echo "  TOOLS were installed in ${bin_nek_tools}                    ";\
	echo "##############################################################";

clean:
	@for i in $(MODULES) ; do         \
      if test -d $${i} ;  then       \
        cd $${i} ;                   \
        ${MAKE} clean ;              \
        cd .. ;                      \
      fi ;                           \
    done ;                           \
	 rm -f $(bin_nek_tools)/nekk      \
			 $(bin_nek_tools)/neks      \
 			 $(bin_nek_tools)/prex      \
 			 $(bin_nek_tools)/nekmerge  \
 			 $(bin_nek_tools)/reatore2  \
 			 $(bin_nek_tools)/genbox    \
  			 $(bin_nek_tools)/genmap    \
			 $(bin_nek_tools)/n2to3;    \
