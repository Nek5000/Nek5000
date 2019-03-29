Reads a CGNS file and generates a .re2 file

* Dependencies: CGNS library 
* Supported element types are HEXA8, HEXA20 and HEXA27
* Boundary IDs are stored in the 5th argument of the fluid bc array in re2
* The real BC's have to be specified in the .usr file using the boundaryID

Current Limitations
-------------------
* Single zone 3D unstructured grid CGNS file using ADF
* Periodic boundary faces are aligned with global coordinate system and match 1:1
