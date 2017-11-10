# Release v17.0-rc1

## Major Features and Improvements

* New paramater file `.par` (replacing .rea)
* Added `OIFS` for `moving mesh`
* Added `Moving mesh` for `PN/PN`
* Improved stability for varying visosity and `PN/PN`
* Added support for mixed `Helmholtz/CVODE` solves
* New fast `AMG setup` tool based on HYPRE
* New `EXODUSII` mesh converter
* New interface to `libxsmm` (fast MATMUL library).
* Extended `lowMach` solver for time varying thermodynamic pressure
* Added DG for scalars
* Added support for `implicit none` in .usr file
* Reduced solver initilization time (parallel binary reader for .re2 and .ma2)
* Restart from arbitrary `fld-file` (multiple files not supported) using interpolation
* Optional new user friendly `SIZE` file format (see SIZE.template)
* Refactored `NEKNEK`
* Added high-pass filter relaxation (alternative to explicit filter)
* Added parallel readers for .re2 and AMG .dat files
* Introduced new binary map (.ma2) format

## Backwards-Incompatible Changes 

* Replaced usr interpolation wrapper `intpts()` by `intp()` with a different interface
* Replaced `g2gi()` by new generic fld reader `gfldr()`
* Moved `makenek` from `core` to `bin` folder
* Removed `MOAB` support 
* Replaced `hpts.in/hpts.out` by `<casename>.his` 
* Eliminated PPLIST symbol `MPIIO` as it is enabled by default now (only active if p65=1 or nfiler=1)
* Eliminated PPLIST symbol `AMG_DUMP` as we dump the files automatically if they don't exist
* Eliminated PPLIST symbol `AMG` as it is a runtime parameter now (rea:p40 or par:preconditoner=semg_amg in PRESSURE section) 
* Changed various key/values for `.par`
* Changed meaning of param(26) to be the targetCFL instead of OIFS substeps

## Known Bugs 

* NONE 

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@ggiannako, @nicooff, @kmittal2, @cliosaglietti, @EvelynOtero, @mattiabr, @maxhutch, @hackljf, @negips 


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
