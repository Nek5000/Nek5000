Reads a GMSH file and generates a .re2 file

* Supported element types are QUAD8, QUAD9, HEXA20 and HEXA27
* Boundary IDs are stored in the 5th argument of the fluid bc array in re2
* The real BC's have to be specified in the .usr file using the boundaryID

Limitations
-------------------
* Unstructured grid GMSH file format 2.x
* A periodic boundary face pair has to match 1:1 
