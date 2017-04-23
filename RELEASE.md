# Release 17.?.?

## Major Features and Improvements

* New paramater file `.par`
* Added `OIFS` for `moving mesh`
* Added `Moving mesh` for `PN/PN`
* Added support for mixed `Helmholtz/CVODE` solves
* New `AMG setup` tool based on HYPRE
* New `EXODUSII` mesh converter
* New interface to `libxsmm` (fast MATMUL library).
* Extended `lowMach` solver for time varying thermodynamic pressure
* Added DG for scalars
* Added support for `implicit none` in .usr file
* Read `.re2` in parallel
* Restart from arbitrary `fld-file` (multiple files not supported) using interpolation
* Optional new user friendly `SIZE` file format (see SIZE.template)
* Refactored `makenek` for faster builds

## Backwards-Incompatible Changes 

* When p20>0 use it as solver tolerance for temperature instead of p22 (only for Helmholtz) 
* Replaced interpolation wrapper `intpts()` by `intp()` with a different interface
* Replaced `g2gi()` by new generic fld reader `gfldr.h`
* Moved `makenek` from core to bin folder
* Removed `MOAB` support 
* Replaced `hpts.in` by `his.in`  and `hpts.out` by `<casename>.his` 
* Eliminated PPLIST symbol `MPIIO` as it is enabled by default now (only active if p65=1 or nfiler=1)
* Eliminated PPLIST symbol `AMG_DUMP` as we dump the files automatically if needed  
* Eliminated PPLIST symbol `AMG` as it is a runtime parameter now (rea:p40 or par:solver=amg) 

## Bug Fixes and Other Changes

* Use rank id to tag MPI message instead of global element 
* Fix periodic BC bug of `genmap` when using more than 1M elemtents
* Many bug fixes

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@ggiannako, @nicooff, @kmittal2, @cliosaglietti, @EvelynOtero, @adampep, @mattiabr, @maxhutch, @hackljf, @negips 


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
