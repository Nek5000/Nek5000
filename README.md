# Nek5000

| **`Short Tests`** | **`Examples`** |
|-----------------|---------------------|
| [![Build](https://travis-ci.org/Nek5000/Nek5000.svg?branch=master)](https://travis-ci.org/Nek5000/Nek5000) | [![Build Status](https://jenkins-ci.cels.anl.gov/buildStatus/icon?job=Nek5000)](https://jenkins-ci.cels.anl.gov/job/Nek5000/) |

Nek5000 is a fast and scalable open source CFD solver. 
Make sure to read the [release notes](https://github.com/Nek5000/Nek5000/blob/master/RELEASE.md) before using the code.
. 

## Directory Structure

Here's a brief description of each top-level directory:

#### `core`
contains the Nek5000 application sources.

#### `bin`
contains scripts for running nek5000 and manipulating its output.

#### `tools`
contains the sources for the pre- and post-processing tools which are stand-alone.

#### `short-tests`
contains light-weight regression tests for verification.

#### `run`
consistent place for users to place their problem cases.

#### `examples`
reserved for the Nek5000 example problems.

#### `doc`
contains the user documentation in HTML and PDF.

#### `3rd_party`
its purpose it to provide a consistent place for 3rd party code.

## Case Files

#### `SIZE`
contains some hardwired runtime parameters to dimension static arrays

#### `foo.par`
contains runtime parameters 

#### `foo.re2`
contains mesh and boundary data

#### `foo.ma2`
contains partioning data

#### `foo.usr`
contains user specific code to initialize solver, set source terms and boundary conditions or to manipulate solver internals. 

#### `foo.his`
contains probing points

#### `foo.f00000`
contains checkpoint data

#### `foo.nek5000`
contains metadata for VisIt

#### `foo.rea` (legacy)
contains runtime parameters and mesh in ASCII. Replaced by .par and .re2 file

#### `foo.map` (legacy)
contains partioning data in ASCII

**Note:** The old legacy files (.rea & .map) are recommended for debugging purposes only. 

## Scripts

Let's walk us through some useful batch scripts:

* `makenek <case>` compiles your case
* `nek/nekb <case>` runs a serial job in foreground or background
* `nekmpi/nekbmpi <case> <number of ranks>` runs a parallel job
* `neknek <case1> <cas2> <ranks 1> <ranks 2>` runs Nek5000 with two overlapping component grids
* `visnek <case>` creates metadata file required by [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
* `mvn <old name> <new name>` renames all case files
* `cpn <old name> <new name>` copies all case files

## Running your very first simulation

Hold your horses, this needs less than 5 min.

```
cd ~
tar -xvzf Nek5000_17.0.0.tar.gz
export PATH=$HOME/Nek5000/bin:$PATH
cd ~/Nek5000/tools; ./maketools genmap
cd ~/Nek5000/run
cp -r ../examples/turbChannel .
cd turbChannel
genmap                    # run partioner, on input type ethier
makenek turbChannel       # build case, edit script to change settings
nekbmpi turbChannel 2     # run Nek5000 on 2 ranks in the background
echo -2 >ioinfo           # stop Nek5000 run and dump checkpoint
visnek turbChannel; visit -o turbChannel.nek5000 # requires a VisIt installation

```

## Meshing

Nek5000 is mainly a solver. However, simple box type meshes can be generated with the `genbox` tool. For more complex meshes please consider using `PRENEK` and the meshing tools `nekmerge` and `n2to3`. We provide mesh converters like `exo2nek` and `msh2nek` which are quite handy if you want to use your favorite mesh generator. Also check our [Bazaar](https://github.com/Nek5000/NekBazaar) for 3rd party meshing tools. 

## Visualization

Nek5000 output (fld) files can be read by [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) or [ParaView](https://www.paraview.org/). There is also an build-in X-Window based postprocessor called `POSTNEK` located in tools.

## Documentation

Visit our online [User's Guide](http://Nek5000.github.io/NekDoc/) which is also available in [PDF](http://nek5000.github.io/NekDoc/Nek_users.pdf).

## Troubleshooting

If you run into problems compiling, installing, or running Nek5000, first check the User's Guide. If you are not able to find a solution to your problem there, please send a message to the User's Group [mailing list](https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users).

## Reporting Bugs
Nek5000 is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/Nek5000/issues) feature on GitHub. However, GitHub Issues should not be used for common troubleshooting purposes. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group mailing list. If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code, you may create an Issue on GitHub.

## Contributing

Our project is hosted on [GitHub](https://github.com/Nek5000). If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

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
5. Set upstream for local master branch `git branch --set-upstream master remotes/origin/master`
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
2. Implement your code changes. To reset your branch and discard any changes run `git reset --hard origin/master`. To revert a set of files run `git checkout file1 file2 ...`
3. Commit your changes to your local repo using e.g. `git commit file1 file2 ...`. Do this frequently to save your work.
4. Periodically, changes made in our master should be pulled back into your local branch by `git pull -r`. This ensures that we do not end up in integration hell that will happen when many feature branches need to be combined at once.
5. If there are no merge conflicts, go to the next step. In case of conflicts edit the unmerged files in question. Merge conflicts are indicated  by the conflict marker `<<<<<<<` in your file.
6. Assuming you are happy run `nekgit_push`. This will create a pull request on GitHub. You can check with `git diff origin/master` what your push will do. When your pull request was merged, run `git pull` on your local master branch to see your change. You can delete the branch created in step (1) with `nekgit_rm <my branch name>`.
