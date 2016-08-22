# CMT-nek #

CMT-nek is an extension of nek5000 to solve compressible flows via the discontinuous Galerkin method.

## Getting started ##

Quickstart for running, with examples for Ubuntu.

### Dependencies ####

-  MPI implementation

>  `  sudo apt-get install libopenmpi-dev openmpi-bin`

-  fortran compiler

>  `  sudo apt-get install gfortran`

-  tsch for some scripts in tools

>  `  sudo apt-get install tcsh`

### Running an example ###

Steps to run an example, in this case swept_particles
- Copy desired example to scratch space
> cp -r /path/to/cmt-nek/examples/CMT/swept_particles /path/to/scratch

- Copy makenek into the example directory
>  `  cp /path/to/cmt-nek/nek/makenek /path/to/scratch/swept\_particles/`

- Edit makenek to set SOURCE_ROOT as the absolute root path to the root of nek

- execute makenek with the casename, where the casename is the usr filename (box.usr)

>  `  ./makenek box`

- Run the case using one of the scripts from tools, such as nekmpi box 4, to run the box case with 4 processes

>  `  /path/to/cmt-nek/tools/scripts/nekmpi box 4`

### Viewing a result ###

Use the script visnek in tools to create a nek5000 header that visit and paraview can use

>  `  /path/to/cmt-nek/tools/scripts/visnek fldbox`