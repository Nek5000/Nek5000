# Nek5000 
![https://travis-ci.org/Nek5000/nek5000](https://travis-ci.org/Nek5000/nek5000.svg?branch=develop)
[![GPLv3 licensed](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/Nek5000/nek5000/develop/LICENSE)

## Wait, where did everything go?

With the move to git, we have decided to reorganize the sources to improve modularity.
Most of what used to be in `trunk/nek` now lives in `core`, including the all-important `makenek`.
We hope the rest of the directories are self describing: you can find things like `nekmpi` and `nekb` in `bin` or `genbox` and `genmap` in `tools`.
The examples have been moved into a seperate repository, [nek5000_examples](https://github.com/Nek5000/nek5000_examples), to keep this one light-weight. 

## Get Nek5000

You can download the latest release of nek5000 as [a zip](https://github.com/Nek5000/nek5000/archive/master.zip) or [a tarball](https://github.com/Nek5000/nek5000/archive/master.tar.gz).
You can also clone the repository with git:
```
git clone https://github.com/Nek5000/nek5000.git -b master --single-branch
```
or even check it out with svn!
```
svn co https://github.com/Nek5000/nek5000.git/branches/master nek5000
```

## Use Nek5000
nek5000 works the same way it used to: build cases with `core/makenek` and run them with a script, e.g. `bin/nekmpi`.
For more information, see the [user guide](https://nek5000.mcs.anl.gov/documentation/).


