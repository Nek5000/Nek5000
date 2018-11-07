Reads an exodus II (.exo) file and generates a .re2 file containing the grid,
curved sides and sideset IDs. 

   - Uses dynamic memory allocation to reduce the memory footprint  
   - Requires the 3rd-party Exodus library (and NetCDF with the suggested modifications by Exodus) 
   - Supported element types are HEX20 (3D) and QUAD8 (2D)
   - HEX27/QUAD9 elements may also work but this is not guaranteed, since there is no standard for this element type in the Exodus library.
   - 2D meshes must be constructed on the z=0 plane
   - Sideset IDs are stored in the 5th argument of the fluid bc array
   - The "real" BC's have to be specified in the .usr file using the sideset IDs 
   - Periodic BC's are not supported yet
   - Conjugate heat transfer (v and t-mesh) is not supported
