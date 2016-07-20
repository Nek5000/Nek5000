Reads a exodus II (.exo) file and generates a .re2 file containing the grid,
curved sides and siteset ids. 

   - Tested with Cubit, but any other package exporting exodusII should work too
   - Sideset ids are stored in the 5th argument of the fluid bc array
   - The "real" BC's have to be specified in the .usr file using the sideset ids 
   - For 2d meshes use z=0
   - Periodic BC's are not supported
   - Conjugate heat transfer (v and t-mesh) is not supported

## BUILD INSTRUCTIONS
There are a few externally developed third-party libraries that are required to build the converter

- [ ]  Download latest netcfd-c release from http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
- [ ]  Modify the following defines in include/netcdf.h (recommended for use in exodus)
```
#define NC_MAX_DIMS     65536 
#define NC_MAX_VARS     524288
#define NC_MAX_VAR_DIMS 8
```
- [ ]   ```./configure --prefix=$HOME/lib CC=gcc --disable-netcdf-4 --disable-fsync --disable-dap```
- [ ]  ```make && make install ``
- [ ]  Download latest exodus release from https://github.com/gsjaardema/seacas/archive/exodus.zip
- [ ]  Ensure ```NetCDF_DIR:PATH``` in ```cmake-exodus points``` points to the netcdf source
- [ ]   ```mkdir build; cd build ```
- [ ]  ```../cmake-exodus ``` 
- [ ]  ```make && make install``` 
- [ ]  Run ```maketools exo2nek``` (add library location to search path using USR_LFLAGS)

