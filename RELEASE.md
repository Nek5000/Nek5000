# Release 17.0.0

## Major Features and Improvements

* New paramater file `.par`
* Added generic `fld` reader to restart from arbitrary mesh
* Added `OIFS` for `moving mesh`.
* Added `Moving mesh` for `PN/PN`.
* Each scalar field can be solved by `Helmholtz` or `CVODE`.
* Improved `.re2` reader (reading in parallel).
* Added `HPF-RT` (high pass filter relaxation term) for stabilization.
* Added `k-omega` turbulence model.
* New `AMG setup` tool based on HYPRE. 
* New `EXODUSII` mesh converter.
* New interface to `libxsmm` (fast MATMUL library).
* Extended `lowMach` solver for time varying thermodynamic pressure.
* Added DG for scalars
* Added support for `implicit none`

## Breaking Changes to the API

* Optional `intp.h` module replaced old interpolation routines `intpts()`
* New generic fld reader `gfldr.h` replaced grid-to-grid interpolation `g2gi()`
* Moved `makenek` to bin folder
* New `SIZE` file with explicit declarations and improved readability

## Bug Fixes and Other Changes

* Use rank id to tag MPI message instead of global element 
* Fix periodic BC bug of `genmap`.
* Many bug fixes

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@ggiannako, @nicooff, @kmittal2, @cliosaglietti, @EvelynOtero, @adampep, @mattiabr, @maxhutch, @hackljf, @negips 


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.