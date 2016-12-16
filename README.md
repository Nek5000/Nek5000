# Nek5000 
![https://travis-ci.org/Nek5000/Nek5000](https://travis-ci.org/Nek5000/Nek5000.svg?branch=develop)
[![GPLv3 licensed](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/Nek5000/nek5000/develop/LICENSE)

Nek5000 is an open source, highly scalable and portable spectral element code designed to simulate:

* unsteady Stokes
* unsteady incompressible Navier-Stokes
* low Mach-number flows
* heat transfer and species transport
* incompressible magnetohydrodynamics (MHD)

Written in Fortran 77 and C, it relies on MPI for parallelization and on a
subset of LAPACK routines for eigenvalue computations, depending on the
particular solver employed.  Nek5000 output formats can be read by the parallel
visualization package VisIt developed at LLNL/LBNL. 


## Features

* High-order spatial discretization using spectral elements
* High-order semi-implicit timestepping
* Incompressible + low Mach number (variable density) flows
* Low memory footprint and scalable memory design
* Scales to over a million processes
* Efficient preconditioners 
* Highly optimized computational kernels 
* High performance parallel I/O
* ALE / moving meshes and free surface flow
* Accurate Lagrangian particle tracking
* Conjugate fluid-solid heat transfer
* Built-in profiling analysis
* Interface to VisIt for parallel data analysis and visualization


## Download

You can download the latest release of nek5000 [here](https://github.com/Nek5000/nek5000/archive/master.tar.gz).


## Getting Started

0. Unpack the tarball e.g. to `~/src/Nek5000`
1. Add `Nek5000/bin` to your search path
2. Copy `Nek5000/short-tests/eddy` to ~/Nek5000/eddy
3. Go to `~/Nek5000/eddy`
4. Copy `~/src/Nek5000/core/makenek` into current directory
5. Compile code  `./makenek eddy` (edit makenek to change build settings)
6. Run your case with `nekmpi eddy`

**Note:** Here you'll find more [examples](https://github.com/Nek5000/NekExamples)

## Documentation

Visit our [website](https://nek5000.mcs.anl.gov/)

## Code Structure

Here's a brief description of each top-level directory:

####`core`
contains the majority of the Nek5000 application sources.

####`jl`
contains gather/scatter communication, interpolation, and preconditioners written in highly general C code.
`jl` used to live in `nek5_svn/trunk/nek`, but is being promoted to the top level to emphasize its library-like relationship to the rest of the source.
In fact, `jl` has been extended externally in [gslib](https://github.com/gslib/gslib), which is used in other projects.

####`bin`
contains scripts for running nek5000 and manipulating its output.

#### `tools`
contains the sources for the pre- and post-processing tools which are stand-alone fortran programs.

#### `short-tests` 
contains light-weight regression tests for validation.  

#### `3rd_party`
contains nothing. Its purpose it to provide a consistent place for 3rd part plugin/toolbox developers to place their code.

## Contributing

First off, thanks for taking the time to contribute! Our project is hosed in [GibHub](https://github.com/Nek5000/Nek5000). To clone our repository just run `git clone https://github.com/Nek5000/nek5000.git`. 

Please branch off and open pull requests to the `develop` branch.
The `master` branch is reserved for releases.

### Basic Workflow
1. create a branch hosting your changes with `nekgit my123 develop`
2. implement your changes
3. commit the changes to your local repo using `git commit`
4. bring in the latest changes by `nekgit_pull` and resolve potential conflicts
5. run `nekgit_push` to create a request on GitHub to merge your changes 

**Note:** A branch should include a consistent and atomic change. Do not mix unrelated changes into a single branch. You can work on multiple in parallel. Just switch between them using `git checkout <branch name>`. Also for a bigger/longer change you may want to interate between step 3 and 4. 
