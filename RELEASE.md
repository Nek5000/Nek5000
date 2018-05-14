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

## Backwards-Incompatible Changes 

* Eliminated PPLIST symbol `NEKNEK` (now a runtime parameter)
* New `map/ma2` format (old versions are NOT supported)
* Changed interface of interpolation wrapper, see `intp.f`

## Known Bugs 

[474](https://github.com/Nek5000/Nek5000/issues/474)
[407](https://github.com/Nek5000/Nek5000/issues/407)
[65](https://github.com/Nek5000/Nek5000/issues/65)

## Thanks to our Contributors
This release contains contributions from the Nek5000 core developers, as well as:

@ggiannako, @kmittal2, @RonRahaman, @kentO


We are also grateful to all who filed issues or helped resolve them, asked and answered questions, and were part of inspiring discussions.
