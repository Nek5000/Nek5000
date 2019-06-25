# Release v19.0-rc2

## What is new? 

* Uncoupled multisession (neknek) simulations
* Gather scatter operations across sessions
* Gather scatter options across gtp-planes
* par file support for postnek
* mkSIZE to automatically create SIZE file
* Object and boundary handling
* RANS k-Omega and k-Omega-SST (experimental) 
* Online mesh-smoother (experimental)
* ElapsedTime option for writeControl (in par)
* Print runtime-statistics every 500 steps
* Support for GNU > 8.x compilers
* Support for Cray compilers
* Support for ARM compilers
* Add AVM regularization for scalars (experimental)
* FEM_AMG precoditioner (experimental) p40=3
* SEMG_AMG_HYPRE precoditioner (experimental) p40=2
* CHT support for generic fld reader
* Overwrite core routines in usr
* Lagrangian phase model - LPM (experimental)
* Distributed gllnid/gllel (experimental) 
* Add parMetis partitioner
* Add parRSB partitioner
* Add CGNS mesh converter
* Support p0th with Helmholtz 
* Update to GSLIB v1.0.5 and HYPRE v2.15.1
* Various bug fixes

## What you may have to change to be compatible 

* Remove PPLIST symbol NEKNEK (not required anymore)
* Use valint instead of ubc in userbc for neknek
* Remove multimesh_create call from usr file (not required anymore)
* Adjust calls to interpolation wrapper according to new interface in interp.f
* Remove common block CTORQ from usr (now part of OBJDATA included in TOTAL)
* Use nekamg_setup tool instead of amg_hypre (required for semg_amg preconditioner) 
* Your parameters to the reserved user space param(170) - param(200) 
* Set lelr in SIZE for restart using muliple files (check value in file header) 
* Use planar_avg() instead of planar_average_z etc. 
* Rename AMG input files (example: amg_Aff.dat -> ethier.amgAff.dat) 

## Known Bugs 

[635](https://github.com/Nek5000/Nek5000/issues/635)
[634](https://github.com/Nek5000/Nek5000/issues/634)
[628](https://github.com/Nek5000/Nek5000/issues/628)
[627](https://github.com/Nek5000/Nek5000/issues/627)
[562](https://github.com/Nek5000/Nek5000/issues/562)
[65](https://github.com/Nek5000/Nek5000/issues/65)

## Thanks to our Contributors

We are grateful to all who added new features, filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
