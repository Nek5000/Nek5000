============================
SERIAL AMG SETUP USING HYPRE
============================

Serial version of the AMG setup for Nek5000 using Hypre, a linear algebra library (http://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods). 
Given the dump files 'amgdmp_i.dat', 'amgdmp_j.dat' and 'amgdmp_p.dat' produced by Nek5000 when 'IFAMG_DUMP=true', the code performs the AMG setup and produces the output files 'amg_Aff.dat', 'amg_AfP.dat', 'amg_W.dat' and 'amg.dat' for running Nek5000 with 'IFAMG=true'.

Hypre
-----

The code requires the use of the Hypre library for performing the AMG setup. Make sure that the library is properly built on your machine before performing the setup:

- Downlad version 2.11.1 from http://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/software

- Extract the files and read install instructions from the file hypre-2.11.1/INSTALL

Makefile
--------

In the Makefile:
- Choose your MPI c compiler. Default is mpicc. Even though the code is serial, Hypre requires MPI.

- Set path to Hypre directory: HYPRE_DIR='PathToHypre'/hypre-2.11.1/src/hypre

- Set compiling flags at your convenience.

Compile
-------

Compile the code typing 'make'.

Run setup
---------

- Copy the 'amgdmp_*.dat' files to your working directory.

- Run the setup by typing './amg_hypre'

- You will be prompted for
    * a coarsening method. Enter the value corresponding to your choice. More info can be found in Hypre user's guide.

    * a maximum number of coarsening levels. Use low values only if you intentionally need to limit the number of levels.

    * a smoother tolerance. This value impacts directly the number of Chebyshev iterations that will be performed during the AMG solve. As a rule of thumb, the tolerance should lie in the range 0.1~1. A lower value leads to more Chebyshev iterations (slower but more accurate), while a high value has the opposite effect.
