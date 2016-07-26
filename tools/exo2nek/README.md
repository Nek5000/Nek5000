Reads a exodus II (.exo) file and generates a .re2 file containing the grid,
curved sides and siteset ids. 

   - Requires the 3rd-party Exodus library (and NetCDF with modifications see below) 
   - So far tested with Cubit, but any other package exporting exodusII should work, too
   - Supported are HEX27 (3D) and QUAD9 (2D)
   - For 2d meshes use z=0
   - Sideset ids are stored in the 5th argument of the fluid bc array
   - The "real" BC's have to be specified in the .usr file using the sideset ids 
   - Periodic BC's are not supported
   - Conjugate heat transfer (v and t-mesh) is not supported

## BUILD INSTRUCTIONS
1.  Download latest netcfd-c release from http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
2.  Modify the following defines in include/netcdf.h (recommended for use in exodus)
    #define NC_MAX_DIMS     65536 
    #define NC_MAX_VARS     524288
    #define NC_MAX_VAR_DIMS 8 
3.  ./configure --prefix=$HOME/lib CC=gcc --disable-netcdf-4 --disable-fsync --disable-dap
4.  make && make install


1. Download latest exodus release from https://github.com/gsjaardema/seacas/archive/exodus.zip
3. Ensure NetCDF_DIR:PATH in cmake-exodus points to the netcdf source
4. mkdir build; cd build;
5. ../cmake-exodus  
6. make && make install
7. Copy exodusII.inc from the exodus source the convert source path
8. Run maketools 

