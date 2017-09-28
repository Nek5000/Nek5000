Reads a exodus II (.exo) file and generates a .re2 file containing the grid,
curved sides and sideset ids. 

   - Requires the 3rd-party Exodus library (and NetCDF with modifications see below) 
   - So far tested with Cubit, but any other package exporting exodusII should work, too
   - Supported are HEX27 (3D) and QUAD9 (2D)
   - For 2d meshes use z=0
   - Sideset ids are stored in the 5th argument of the fluid bc array
   - The "real" BC's have to be specified in the .usr file using the sideset ids 
   - Periodic BC's are not supported yet
   - Conjugate heat transfer (v and t-mesh) is not supported
