# Nek5000 
![https://travis-ci.org/Nek5000/Nek5000](https://travis-ci.org/Nek5000/Nek5000.svg?branch=develop)
[![GPLv3 licensed](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/Nek5000/nek5000/develop/LICENSE)

Nek5000 is an open source, highly scalable and portable spectral element code designed to simulate:

* unsteady Stokes
* unsteady incompressible Navier-Stokes
* low Mach-number flows
* heat transfer and species transport
* incompressible magnetohydrodynamics (MHD)

Written in Fortran 77 and C, it relies on MPI for parallelization. Nek5000 output formats can be read by the parallel visualization package VisIt developed at LLNL/LBNL. 


## Features

* High-order spatial discretization using spectral elements
* High-order semi-implicit timestepping
* Incompressible + low Mach number flows
* RANS and LES turbulence models (experimental)
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

You can download the latest release of Nek5000 [here](https://github.com/Nek5000/nek5000/archive/master.tar.gz).


## Getting Started

1. Unpack the tarball to `~/Nek5000`
2. Add `~/Nek5000/bin` to your shell search path
3. Copy `~/Nek5000/short-tests/eddy` to `~/nekcases/eddy`
4. Copy `~/Nek5000/core/makenek ` to `~/nekcases/eddy`
5. Go to `~/nekcases/eddy` and run `./makenek eddy` (see makenek for build options)
5. You can run the case using two processes with `nekmpi eddy 2`

**Note:** Here you'll find more [examples](https://github.com/Nek5000/NekExamples)

## Documentation

Visit our [website](https://nek5000.mcs.anl.gov/documentation)

## Troubleshooting

If you run into problems compiling, installing, or running Nek5000, first check the [User's Guide](http://nek5000.github.io/NekDoc/Nek_users.pdf). If you are not able to find a solution to your problem there, please send a message to the User's Group [mailing list](https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users).

## Reporting Bugs
Nek5000 is hosted on GitHub and all bugs are reported and tracked through the Issues feature on GitHub. However, GitHub Issues should not be used for common troubleshooting purposes. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group mailing list. If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code, you may create an Issue on GitHub.

## Contributing

First off, thanks for taking the time to contribute! Our project is hosted on [GibHub](https://github.com/Nek5000/Nek5000). The main repository will always hold two evergreen branches:

* `develop`
* `master`

The main branch should be considered `develop` and will be the main branch where the source code of `HEAD` always reflects a state with the latest delivered development changes for the next release. As a developer, you will you typically be branching and merging from `develop`.

Consider `master` to always represent the latest code deployed to production. During day to day development, the `master` branch will not be interacted with. When the source code in the `develop` branch is stable and has been deployed, all of the changes will be merged into `master` and tagged with a release number. 

### One Time Setup
1. Sign up on [GibHub](https://github.com/)
2. Fork our [project](https://github.com/Nek5000/Nek5000) on GitHub
3. Download fork with `git clone -o myfork https://github.com/<username>/Nek5000.git ~/Nek5000`
4. Add our repo `git remote add origin https://github.com/Nek5000/Nek5000.git`
5. Download our repo `git fetch origin`
6. Set upstream for local develop branch `git branch --set-upstream-to origin/develop develop`
7. Run `~/Nek5000/bin/git-hub setup —u <your username on GitHub> --global`
8. Add this to your [hub] section in `~/.gitconfig`:

```
[hub]
        ...
        upstream = Nek5000/Nek5000
        forkremote = myfork 
``` 

### How It Works
1. Create a branch hosting your changes with `nekgit_co <my branch name> develop`. The core idea is that all development should take place in a _dedicated_ branch instead of the local development branch.
2. Implement your changes. Make sure your change is atomic and consistent. You can work on multiple branches simultaneously. Just do a `git checkout <your branch name>` to change the branch. Note, this will update the files in your working directory (~/Nek5000). To compare your branch with our develop repo use `git diff origin/develop`.
3. Commit the changes to your local repo using `git commit -a -m 'a descriptive comment'`. Do this frequently to save your work (otherwise you cannot switch branches). 
4. Periodically, changes made in our Nek5000 repo should be pulled back into your branch by `git pull`.
5. If there are no merge conflicts, go to the next step. In case of conflicts edit the unmerged files in question. Merge conflicts are indicated  by the conflict marker `<<<<<<<` in your file. If you are done with all files, run `git add .` and do a `git commit` to indicate that all conflicts have been resolved.  
6. Assuming you are happy with your change, run `nekgit_push` to create a request on GitHub to merge your changes. Now you should be able to see your pull request on GitHub. The core developers will review your change and discuss potential modifications. We cannot consider your merge request if it is outdated or does not pass the regression tests. Please include a short-test in case of a new feature. When your pull request was merged or closed, you can delete your branch (created in step 1) with `nekgit_rm <my branch name>`.
7. You may want to set your working directory to the latest develop branch. To do this just run `git checkout develop; git pull`. After your pull reqest was merged, you have to update your local develop branch again (git pull) to see your change. 


## Code Structure

Here's a brief description of each top-level directory:

####`core`
contains the majority of the Nek5000 application sources.

####`jl`
contains gather/scatter communication ([gslib](https://github.com/gslib/gslib)), interpolation, and preconditioners written in highly general C code.

####`bin`
contains scripts for running nek5000 and manipulating its output.

#### `tools`
contains the sources for the pre- and post-processing tools which are stand-alone.

#### `short-tests` 
contains light-weight regression tests for validation.  

#### `3rd_party`
contains nothing. Its purpose it to provide a consistent place for 3rd part plugin/toolbox developers to place their code.


