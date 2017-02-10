#AMG HYPRE

Read the AMG data files and produce solver input files.  

#INSTALLATION

To compile the code you need to link against the [HYPRE](http://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) library. Make sure to build HYPRE with `--without-MPI` as the tool works in serial only.

Edit `maketools`

* Add the HYPRE include path to `CC` e.g. `CC=gcc -I/src/hypre/include`
* Use `USR_LFLAGS` to link against HYPRE e.g. `USR_LFLAGS=-L$HOME/src/hypre/lib -lHYPRE`

Then, just run `maketools amg_hypre`


#Workflow

1. Compile and run Nek5000 with the `AMG` preprocessor symbol. This will create the AMG data files `amgdmp_*.dat`
2. Run `amg_hypre` in the same directory to procude the solver input files `amg.dat` and `amg_*.dat` 
3. Finally, run Nek5000 again. 

