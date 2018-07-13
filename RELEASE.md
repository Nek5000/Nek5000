# Release v18.0-3

## Major Features and Improvements

* Introduced uncoupled multisession simulations
* Added gather scatter operations across sessions
* Added gather scatter options across gtp-planes
* Added par file support for postnek
* Added mkSIZE script to automatically create SIZE file
* Refactored object and boundary handling
* Added RANS models (experimental) 
* Added mesh-smoother (experimental)
* Added elapsedTime option for writeControl (in par)
* Print runtime-statistics every 100 steps
* Added support for GNU 8.x
* Added FEM-AMG precoditioner (experimental)
* Added support for CHT to generic fld reader

## Backwards-Incompatible Changes 

* Remove PPLIST symbol `NEKNEK` (not required anymore)
* Use valint instead of ubc in `userbc` for neknek
* Remove multimesh_create call from `usr` file (not required anymore)
* Change to modified interpolation routines (see interp.f)
* Remove CB CTORQ from usr files (now part of OBJDATA included in TOTAL)
* Use amg_setup instead of amg_hypre (tool was renamed) 

## Bug Fixes
[511](https://github.com/Nek5000/Nek5000/issues/511),
[497](https://github.com/Nek5000/Nek5000/issues/497),
[470](https://github.com/Nek5000/Nek5000/issues/470),
[467](https://github.com/Nek5000/Nek5000/issues/467),
[463](https://github.com/Nek5000/Nek5000/issues/463)

## Known Bugs 

[507](https://github.com/Nek5000/Nek5000/issues/507),
[474](https://github.com/Nek5000/Nek5000/issues/474),
[407](https://github.com/Nek5000/Nek5000/issues/407),
[65](https://github.com/Nek5000/Nek5000/issues/65)

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@perrunchin, @ggiannako, @kmittal2, @kentO


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
