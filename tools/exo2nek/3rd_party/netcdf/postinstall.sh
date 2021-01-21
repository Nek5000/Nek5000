#!/bin/bash
# Script to automatically build, install netcdf-fortran.
# Very rought draft.

DOCMAKE=""
DOAUTOTOOL=""
DOACTION=""
TARGVERSION="v4.4.4"

MARG=""
if [ $# -lt 1 ]; then
    echo ""
    echo "WARNING! THIS SCRIPT IS NOT MEANT TO BE RUN MANUALLY."
    echo "WARNING! THIS SCRIPT IS NOT MEANT TO BE RUN MANUALLY."
    echo ""
    exit 1
fi


##
# Check for 'git', exit if it's not found.
##
hash git 2>/dev/null1
if [ $? -eq 1 ]; then
    echo "ERROR: 'git' is required to install netcdf-fortran automatically"
    echo "through this method. Please reinstall git and try again."
    exit 1
fi

while getopts "a:t:" o; do
    case "${o}" in
        t)
            MARG=${OPTARG}
            ;;
        a)
            DOACTION=${OPTARG}
            ;;
        *)
            echo "Specify type with -t. Types are autotools or cmake."
            exit 1
            ;;
    esac
done


###
# Make sure we're specifying an allowed
# build system.
###
case ${MARG} in
    cmake)
        DOCMAKE="TRUE"
        ;;
    autotools)
        DOAUTOTOOL="TRUE"
        ;;
    *)
        echo "Illegal type. Types are autotools or cmake."
        exit 1
        ;;
esac

###
# Make sure we're performing a valid action.
###
case ${DOACTION} in
    build)
        ;;
    install)
        ;;
    *)
        echo "Illegal action."
        exit 1
        ;;
esac


###
# Fetch netcdf-fortran from git if need be.
###

if [ ! -d "netcdf-fortran" ]; then
    git clone http://github.com/unidata/netcdf-fortran

fi

cd netcdf-fortran
git checkout $TARGVERSION

###
# Invoke cmake to build netcdf-fortran
###

if [ "x$DOCMAKE" = "xTRUE" ]; then

    mkdir -p build
    cd build

    if [ "x$DOACTION" = "xbuild" ]; then
        cmake .. -DCMAKE_PREFIX_PATH=@CMAKE_INSTALL_PREFIX@ -DCMAKE_INSTALL_PREFIX=@CMAKE_INSTALL_PREFIX@ -DBUILD_SHARED_LIBS=@BUILD_SHARED_LIBS@ -DTEST_PARALLEL=@ENABLE_PARALLEL@ &&
        make && make test
    fi

    if [ "x$DOACTION" = "xinstall" ]; then
        make install
    fi

fi

if [ "x$DOAUTOTOOL" = "xTRUE" ]; then

    if [ "x$DOACTION" = "xbuild" ]; then
        if [ ! -f "configure" ]; then
            autoreconf -if
        fi

        STATIC_BUILD="--disable-static"
        if [ "xyes" = "xyes" ]; then
            STATIC_BUILD="--enable-static"
        fi

        SHARED_BUILD="--disable-shared"
        if [ "xyes" = "xyes" ]; then
            SHARED_BUILD="--enable-shared"
        fi

        LIBS="-lm " CFLAGS="-I/nfs2/yhaomin2007/new_exo2nek_pull_request/Nek5000-1/tools/exo2nek/3rd_party/netcdf/install/include" LDFLAGS="-L/nfs2/yhaomin2007/new_exo2nek_pull_request/Nek5000-1/tools/exo2nek/3rd_party/netcdf/install/lib" ./configure --prefix=/nfs2/yhaomin2007/new_exo2nek_pull_request/Nek5000-1/tools/exo2nek/3rd_party/netcdf/install $STATIC_BUILD $SHARED_BUILD
        LIBS="-lm " make
        LIBS="-lm " make check
    fi

    if [ "x$DOACTION" = "xinstall" ]; then
        make install
    fi

fi
