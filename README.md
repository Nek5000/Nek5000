# CMT-nek #

CMT-nek is ...


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

- Copy makenek into the example directory

>  `  cp makenek ../cmt\_examples/swept\_particles/`

- Edit makenek to set SOURCE_ROOT as the absolute root path to the root of nek

- execute makenek with the casename, where hte casename is the usr filename (box.usr)

>  `  ./makenek box`

- Run the case using one of the scripts from tools, such as nekmpi box 4, to run the box case with 4 processes

>  `  ../../tools/scripts/nekmpi box 4`

### Viewing a result ###

Use the script visnek in tools to create a nek5000 header that visit and paraview can use

>  `  ../../tools/scripts/visnek fldbox`


