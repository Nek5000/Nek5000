# Release v18.0-4

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
* Print runtime-statistics every 100 steps
* Support for GNU 8.x compilers
* Support for Cray compilers
* Support for ARM compilers
* FEM_AMG precoditioner (experimental) p40=3
* SEMG_AMG_HYPRE precoditioner (experimental) p40=2
* CHT support for generic fld reader
* Overwrite core routines in usr
* Fully scaleable memory footprint PPLIST=DPROCMAP (experimental)
* Online genmap PPLIST=PARRSB (experimental) 
* Various bug fixes

## What you may have to change to be compatible 

* Remove PPLIST symbol NEKNEK (not required anymore)
* Use valint instead of ubc in userbc for neknek
* Remove multimesh_create call from usr file (not required anymore)
* Adjust calls to interpolation wrapper according to new interface in interp.f
* Remove common block CTORQ from usr (now part of OBJDATA included in TOTAL)
* Use amg_setup tool instead of amg_hypre (required for semg_amg preconditioner) 
* Your parameters to the reserved user space param(170) - param(200) 
* Set lelr in SIZE to number of elements per restart file (for muli-file I/O only)
* Use planar_avg() instead of planar_average_z etc. 

## Known Bugs 

[507](https://github.com/Nek5000/Nek5000/issues/507),
[474](https://github.com/Nek5000/Nek5000/issues/474),
[407](https://github.com/Nek5000/Nek5000/issues/407),
[65](https://github.com/Nek5000/Nek5000/issues/65)

## Thanks to our Contributors

We are grateful to all who added new features, filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
