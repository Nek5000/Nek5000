Reads an exodus II (.exo) file and generates a .re2 file containing the grid,
curved sides and sideset IDs. 

* Requires the 3rd-party Exodus library and NetCDF 
* Supported element types are HEX20 (3D) and QUAD8 (2D)
* 2D meshes must be constructed on the z=0 plane
* Sideset IDs are stored in the 5th argument of the fluid bc array in re2
* The real BC's have to be specified in the .usr file using the sideset IDs

Current Limitations
------------------- 
* Periodic BC's are not supported
* Conjugate heat transfer (v and t-mesh) is not supported
