#Build Instructions

Run `./install` (you may need to modify some settings to your environment)

#How to use CVODE

Edit `makenek` and add

```
USR_LFLAGS+=" -L$HOME/Nek5000/3rd_party/cvode/lib"
USR_LFLAGS+=" -lsundials_fcvode -lsundials_cvode"
USR_LFLAGS+=" -lsundials_fnvecparallel -lsundials_nvecparallel" 
```

an compile Nek5000. Then set `solver=cvode` in your `.par` for every `SCALAR` you want to solve with CVODE.