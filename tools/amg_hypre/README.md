#AMG HYPRE

Read the AMG data files and produce solver input files.  

#Build Instructions

* Run `hypre/install`
* Edit `maketools` and add include path to `CC` e.g. `CC=gcc -I./hypre/include`
* Use `USR_LFLAGS` to link against HYPRE e.g. `USR_LFLAGS=-L./hypre/lib -lHYPRE`
* Finally run `maketools amg_hypre`


#Workflow

1. Run Nek5000 using the AMG pressure solver (this is controlled by a runtime parameters) to obtain the required setup files
2. Run AMG setup tool `amg_hypre` to procude the AMG solver input files `amg.dat` and `amg_*.dat` 
3. Finally, just run Nek5000 again (make sure you don't delete these files accidentally) 
