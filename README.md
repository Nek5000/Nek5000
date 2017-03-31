# Nek5000

| **`Short Tests`** | **`Examples`** |
|-----------------|---------------------|
| [![Build](https://travis-ci.org/Nek5000/Nek5000.svg?branch=develop)](https://travis-ci.org/Nek5000/Nek5000) | [![Build Status](https://jenkins-ci.cels.anl.gov/buildStatus/icon?job=Nek5000)](https://jenkins-ci.cels.anl.gov/job/Nek5000/) |

Nek5000 is an open source, fast and scalable spectral element CFD code designed to simulate unsteady incompressible and low Mach-number Navier-Stokes flows.

## Highlights

* Runs on all POSIX compliant operating systems
* Written in Fortran77 and C
* Pure MPI for parallelization
* Proven scalability to over a million ranks
* Easy-to-build with minimal dependencies
* High-order conformal curved quadrilateral/hexahedral meshes
* 2nd/3rd order adaptive semi-implicit timestepping
* Efficient multigrid preconditioners
* Parallel I/O
* Moving mesh and free surface flow
* Lagrangian particle tracking
* Magnetohydrodynamics (MHD)
* Conjugate fluid-solid heat transfer
* uRANS and LES turbulence models
* Meshing tools and converters
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit) & [Paraview](http://www.paraview.org/) support for data analysis and visualization


## Download

For a typical user we recommend to download the [latest release](https://github.com/Nek5000/nek5000/archive/tbd.tar.gz) once it available. All developers please read the `Contributing` section below.

## Getting Started

Hold your horses in less than 5min you have performed your first simulation

```
cd ~
tar -xvzf Nek5000.tar.gz
export PATH=~/Nek5000/bin:$PATH
cd ~/Nek5000/tools; ./maketools genmap
cd ~/Nek5000/run; cp -r ~/Nek5000/short-tests/eddy .
cp /Nek5000/core/makenek .
./makenek eddy_uv
genmap             # on input type eddy_uv
nekmpi eddy_uv 2   # to run on 2 ranks
```

**Note:** For more information see [here](http://nek5000.github.io/NekDoc/Nek_usersch2.html)

## Documentation

Visit our [website](https://nek5000.mcs.anl.gov/documentation)

## Troubleshooting

If you run into problems compiling, installing, or running Nek5000, first check the [User's Guide](http://nek5000.github.io/NekDoc/Nek_users.pdf). If you are not able to find a solution to your problem there, please send a message to the User's Group [mailing list](https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users).

## Reporting Bugs
Nek5000 is hosted on GitHub and all bugs are reported and tracked through the Issues feature on GitHub. However, GitHub Issues should not be used for common troubleshooting purposes. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group mailing list. If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code, you may create an Issue on GitHub.

## Contributing

Our project is hosted on [GibHub](https://github.com/Nek5000/Nek5000). If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

### How we do it
- Anything in master is always deployable
- Upcoming feature release get their own tags or branch that are branched out of master
- All development happens on the master branch.
- To work on something new, create a short lived local branch off of master
- When you need feedback or help, or your change is ready for merging, open a pull request.

### One-time Setup
1. Fork our [GitHub project](https://github.com/Nek5000/Nek5000)
2. Download fork with `git clone -o myfork https://github.com/<username>/Nek5000.git ~/Nek5000`
3. Add our repo `cd ~/Nek5000; git remote add origin https://github.com/Nek5000/Nek5000.git`
4. Download our repo `git fetch origin`
5. Set upstream for local develop branch `git branch --set-upstream develop remotes/origin/develop`
6. Run `~/Nek5000/bin/git-hub setup —u <your username on GitHub> --global`
7. Add this to your [hub] section in `~/.gitconfig`:

```
[hub]
        ...
        upstream = Nek5000/Nek5000
        forkremote = myfork
```

### Typical Workflow
1. Create a feature branch hosting your change with `nekgit_co <descriptive name>`. Using a dedicated branch for every feature helps you to move between different developments while some are work in progress or under review.
2. Implement your code changes. To reset your branch and discard any changes run `git reset --hard origin/master`. To revert a set of files run `git checkout <file1 file2 ...>`
3. Commit your changes to your local repo using e.g. `git commit -a`. Do this frequently to save your work.
4. Periodically, changes made in our master should be pulled back into your local branch by `git pull -r`. This ensures that we do not end up in integration hell that will happen when many feature branches need to be combined at once.
5. If there are no merge conflicts, go to the next step. In case of conflicts edit the unmerged files in question. Merge conflicts are indicated  by the conflict marker `<<<<<<<` in your file.
6. Check with `git diff origin/master` what your push will do. Assuming you are happy run `nekgit_push`. This will create a pull request on GitHub and set your current branch back to master. When your pull request was merged, make sure you are on your local master branch. Then, delete the branch created in step (1) with `nekgit_rm <my branch name>` and update the current master branch using `git pull`.

## Code Structure

Here's a brief description of each top-level directory:

#### `core`
contains the majority of the Nek5000 application sources.

#### `jl`
contains gather/scatter communication ([gslib](https://github.com/gslib/gslib)), interpolation, and preconditioners written in highly general C code.

#### `bin`
contains scripts for running nek5000 and manipulating its output.

#### `tools`
contains the sources for the pre- and post-processing tools which are stand-alone.

#### `short-tests`
contains light-weight regression tests for validation.

#### `run`
contains nothing. Its purpose it to provide a consistent place for users to place their cases.

#### `3rd_party`
contains nothing. Its purpose it to provide a consistent place for 3rd part developers to place their code.
