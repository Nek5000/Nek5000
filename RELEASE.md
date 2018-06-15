# Release v18.0-1

## Major Features and Improvements

* Reduced memory usage for large element counts
* Improved interpolation performance on low-order geometries
* Introduced uncoupled multisession simulations
* Added gather scatter operations across sessions
* Added gather scatter options across gtp-planes
* Added par file support for postnek
* Added parallel map file reader
* Added mkSIZE script to automatically create SIZE file
* Refactored object and boundary handling
* Added RANS and mesh-smoother as experimental features

## Backwards-Incompatible Changes 

* Eliminated PPLIST symbol `NEKNEK` (not required anymore)
* `map/ma2` format has changed, rerun with shipped genmap
* Modified interpolation routine names and interface
* Remove CB CTORQ from usr files (now part of OBJDATA included in TOTAL)

## Known Bugs 

[474](https://github.com/Nek5000/Nek5000/issues/474)
[407](https://github.com/Nek5000/Nek5000/issues/407)
[65](https://github.com/Nek5000/Nek5000/issues/65)

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@ggiannako, @kmittal2, @RonRahaman, @kentO


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
