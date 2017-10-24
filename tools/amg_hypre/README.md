#AMG HYPRE

Read the AMG data files and produce solver input files.  

#Workflow

1. Run Nek5000 using the SEMG_AMG preconditioner (through .rea or .par) to obtain the required setup files
2. Run AMG setup tool `amg_hypre` to procude the AMG solver input files `amg.dat` and `amg_*.dat` 
3. Finally, just run Nek5000 again (make sure you don't delete these files accidentally) 
