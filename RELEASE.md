# Release 17.?.?

## Major Features and Improvements

* New paramater file `.par`
* Added generic `fld` reader to restart from arbitrary mesh
* Added `OIFS` for `moving mesh`.
* Added `Moving mesh` for `PN/PN`.
* Added support for mixed `Helmholtz/CVODE` solves.
* Improved `.re2` reader (reading in parallel).
* New `AMG setup` tool based on HYPRE. 
* New `EXODUSII` mesh converter.
* New interface to `libxsmm` (fast MATMUL library).
* Extended `lowMach` solver for time varying thermodynamic pressure.
* Added DG for scalars
* Added support for `implicit none`
* New `generic fld` reader allows restarts from an arbitrary mesh

## Backwards-Incompatible Changes 

* Optional `intp.h` module replaced old interpolation routines `intpts()`
* Replcaed `g2gi()` by new generic fld reader `gfldr.h`
* Moved `makenek` to bin folder
* New `SIZE` file required to use `implicit none`
* Eliminated PPLIST symbol `AMG_DUMP` as we dump the files automatically if needed  
* Removed `MOAB` 
* Replaced `hpts.in` & `hpts.out` by `his.in` & `<casename>.his` 


## Bug Fixes and Other Changes

* Use rank id to tag MPI message instead of global element 
* Fix periodic BC bug of `genmap` when using more than 1M elemtents
* Many bug fixes

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@ggiannako, @nicooff, @kmittal2, @cliosaglietti, @EvelynOtero, @adampep, @mattiabr, @maxhutch, @hackljf, @negips 


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
