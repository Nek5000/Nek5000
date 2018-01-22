### makefile automatically created by makenek 01/18/2018 18:50:39 ###
BINNAME=nek5000
LIBNAME=libnek5000.a
CASENAME:=eddy_dir443
CASEDIR:=/home/e/eotero/pfs/Nek5000_acc/Nek5000/short_tests/eddy3d
S:=/home/e/eotero/pfs/Nek5000_acc/Nek5000/core
SCMT:=/home/e/eotero/pfs/Nek5000_acc/Nek5000/core/cmt
J:=/home/e/eotero/pfs/Nek5000_acc/Nek5000/jl
OPT_INCDIR=./
OBJDIR=obj
IFMPI:=false
IFCMT:=false
IFAMG:=false
IFNEKNEK:=false
FC:=pgf90
CC:=pgcc
FFLAGS:=-r8 -Mpreprocess
CFLAGS:=-DNEW_GS_LOOPS
PPPO=
PPS= NOMPIIO UNDERSCORE TIMER GLOBAL_LONG_LONG
G:=-g -mcmodel=medium -acc -Minfo=accel -ta=nvidia
USR:=
USR_LFLAGS:=-g -mcmodel=medium -acc -Minfo=accel -ta=nvidia
IFVISIT=
VISIT_IFLAGS=
VISIT_LFLAGS=

################################################################################

LDFLAGS = $(USR_LFLAGS)

PPS_F = $(patsubst %,$(PPPO)-D%,$(PPS))
PPS_C = $(patsubst %,-D%,$(PPS))

# CORE##########################################################################
CORE = drive.o drive1.o drive2.o \
plan4.o bdry.o coef.o conduct.o connect1.o connect2.o \
dssum.o edgec.o eigsolv.o gauss.o genxyz.o navier1.o makeq.o \
navier0.o navier2.o navier3.o navier4.o prepost.o speclib.o \
map2.o turb.o mvmesh.o ic.o ssolv.o planx.o math.o mxm_wrapper.o \
hmholtz.o gfdm_par.o  gfdm_op.o gfdm_solve.o subs1.o subs2.o \
genbox.o gmres.o hsmg.o convect.o induct.o perturb.o \
navier5.o navier6.o navier7.o navier8.o fast3d.o fasts.o calcz.o \
byte.o chelpers.o byte_mpi.o postpro.o \
cvode_driver.o nek_comm.o \
vprops.o makeq_aux.o \
papi.o nek_in_situ.o \
readat_new.o reader_re2.o finiparser.o iniparser.o dictionary.o acc.o
################################################################################
# MXM 
MXM=mxm_std.o
################################################################################
# BLAS
BLAS=blas.o dsygv.o ssygv.o
# JL Routines ###################################################################
JO  = jl_

JLCORE = $(JO)gs.o $(JO)sort.o $(JO)sarray_transfer.o $(JO)sarray_sort.o \
$(JO)gs_local.o $(JO)crystal.o $(JO)comm.o $(JO)tensor.o $(JO)fail.o \
$(JO)fcrystal.o

# JL INTERPOLATION
JLINTP = $(JO)findpts.o $(JO)findpts_local.o \
$(JO)obbox.o $(JO)poly.o $(JO)lob_bnd.o \
$(JO)findpts_el_3.o $(JO)findpts_el_2.o

# YP NEKNEK SOLVER
ifeq ($(IFNEKNEK),true)
NEKNEKO := multimesh.o
else
NEKNEKO := singlmesh.o
endif

# JL CRS GRID SOLVER
ifeq ($(IFAMG),true)
CGS := $(JO)amg.o $(JO)fcrs.o
else
CGS = $(JO)sparse_cholesky.o $(JO)xxt.o $(JO)fcrs.o
endif

JL := -DPREFIX=jl_

COMM_MPI := comm_mpi.o
ifeq ($(IFMPI),false)
  COMM_MPI := ${COMM_MPI} mpi_dummy.o
endif

# IFCMT append *.o files to CORE
ifeq ($(IFCMT),true)
  CMTCORE := drive1_cmt.o drive2_cmt.o driver3_cmt.o \
   face.o ausm.o  MixtPerf.o intpdiff.o diffusive_cmt.o eqnsolver_cmt.o surface_fluxes.o \
   bc.o inflow_bc.o outflow_bc.o wall_bc.o step.o \
   filters_cmt.o diagnostics.o
else
  CMTCORE := 
endif
# VISIT ########################################################################
ifeq ($(IFVISIT),true)
  VISITNEK_INCLUDES:=-I$S/3rd_party
  VISITO := visit.o 
  LDFLAGS := ${LDFLAGS} $(VISIT_LFLAGS)
endif

ifeq ($(IFMPI),false)
	DUMMY:= $(shell cp $S/mpi_dummy.h $S/mpif.h) 
else
	DUMMY:= $(shell rm -rf $S/mpif.h) 
endif

TMP1 = $(CORE) $(MXM) $(USR) $(COMM_MPI) $(NEKNEKO) $(VISITO) $(BLAS)

ifeq ($(IFCMT),true)
TMP1 += $(CMTCORE)
OPT_INCDIR += -I $(SCMT)
endif

NOBJS_F = $(patsubst %,$(OBJDIR)/%,$(TMP1))
TMP2 = $(JLCORE) $(JLINTP) $(CGS) 
NOBJS_C = $(patsubst %,$(OBJDIR)/%,$(TMP2))
NOBJS  = $(NOBJS_F) $(NOBJS_C)

L0 = $(G) -O0
L2 = $(G) -O2
L3 = $(G) -O3

FL0   = $(L0) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR) 
FL2   = $(L2) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)
FL3   = $(L3) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)

cFL0   = $(L0) $(CFLAGS) $(PPS_C) -I$J
cFL2   = $(L2) $(CFLAGS) $(PPS_C) -I$J
cFL3   = $(L3) $(CFLAGS) $(PPS_C) -I$J

################################################################################
all : nek5000

objdir: 
	@mkdir $(OBJDIR) 2>/dev/null; cat /dev/null 

nek5000:	objdir $(NOBJS)
	$(FC) -c $(FL2) $(CASEDIR)/${CASENAME}.f $(VISITNEK_INCLUDES) -o ${OBJDIR}/${CASENAME}.o 
	$(FC) -o ${BINNAME} ${OBJDIR}/${CASENAME}.o $(NOBJS) $(LDFLAGS)
	@if test -f ${BINNAME}; then \
        echo ""; \
	echo "#############################################################"; \
	echo "#                  Compilation successful!                  #"; \
	echo "#############################################################"; \
        echo ""; \
        echo "$(NEK_WARN)"; \
	else \
	echo -e "\033[1;31;38m" "ERROR: Compilation failed!"; \
	echo -e "\033[0m"; \
	fi
ifeq ($(IFMPI),false) 
	@rm -rf $S/mpif.h
endif

lib:	objdir $(NOBJS)
	$(AR) cru ${LIBNAME} $(NOBJS)
	ranlib ${LIBNAME}

clean:
	rm -rf ${CASENAME}.f obj ${BINNAME}
ifeq ($(IFMPI),false) 
	@rm -rf $S/mpif.h
endif

$(NOBJS_F) : SIZE

# NEK CORE     ##################################################################
$(OBJDIR)/drive.o	:$S/drive.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive1.o	:$S/drive1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive2.o	:$S/drive2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/prepost.o	:$S/prepost.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/postpro.o	:$S/postpro.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gfldr.o	:$S/gfldr.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/connect1.o	:$S/connect1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/connect2.o	:$S/connect2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/edgec.o	:$S/edgec.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/genxyz.o	:$S/genxyz.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/subs1.o	:$S/subs1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/subs2.o	:$S/subs2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/turb.o	:$S/turb.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/map2.o	:$S/map2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mvmesh.o	:$S/mvmesh.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/induct.o	:$S/induct.f;	      		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/convect.o	:$S/convect.f;	      		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/perturb.o	:$S/perturb.f;	      		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/genbox.o	:$S/genbox.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/hsmg.o	:$S/hsmg.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gmres.o	:$S/gmres.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/papi.o	:$S/papi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/bdry.o	:$S/bdry.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ic.o		:$S/ic.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/conduct.o	:$S/conduct.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier0.o	:$S/navier0.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier2.o	:$S/navier2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier3.o	:$S/navier3.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier4.o	:$S/navier4.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier5.o	:$S/navier5.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier6.o	:$S/navier6.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier7.o	:$S/navier7.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier8.o	:$S/navier8.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/speclib.o	:$S/speclib.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/calcz.o	:$S/calcz.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/fast3d.o	:$S/fast3d.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/dssum.o	:$S/dssum.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/eigsolv.o	:$S/eigsolv.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gauss.o	:$S/gauss.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/planx.o	:$S/planx.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ssolv.o	:$S/ssolv.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gfdm_par.o	:$S/gfdm_par.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gfdm_solve.o	:$S/gfdm_solve.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gfdm_op.o	:$S/gfdm_op.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/coef.o	:$S/coef.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/plan4.o	:$S/plan4.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/vprops.o	:$S/vprops.f;	   		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/cvode_driver.o :$S/cvode_driver.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/makeq.o      	:$S/makeq.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/makeq_aux.o  	:$S/makeq_aux.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/hmholtz.o	:$S/hmholtz.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier1.o	:$S/navier1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/fasts.o	:$S/fasts.f;		        $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/comm_mpi.o	:$S/comm_mpi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mpi_dummy.o	:$S/mpi_dummy.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/byte_mpi.o	:$S/byte_mpi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/math.o	:$S/math.f;			$(FC) -c $(FL3) $< -o $@
$(OBJDIR)/multimesh.o	:$S/multimesh.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/singlmesh.o	:$S/singlmesh.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/lb_setqvol.o	:$S/lb_setqvol.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/readat_new.o	:$S/readat_new.f $S/PARDICT;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/reader_re2.o	:$S/reader_re2.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/acc.o	:$S/acc.f;		$(FC) -c $(FL2) $< -o $@

# CMT CORE     ##################################################################
$(OBJDIR)/step.o	:$(SCMT)/step.f;        	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive1_cmt.o	:$(SCMT)/drive1_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive2_cmt.o	:$(SCMT)/drive2_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/driver3_cmt.o	:$(SCMT)/driver3_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gradient_solver.o	:$(SCMT)/gradient_solver.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/eqnsolver_cmt.o	:$(SCMT)/eqnsolver_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ausm.o	:$(SCMT)/ausm.f;   		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/face.o	:$(SCMT)/face.f;  		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/bc.o		:$(SCMT)/bc.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/surface_fluxes.o	:$(SCMT)/surface_fluxes.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/diffusive_cmt.o	:$(SCMT)/diffusive_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/MixtPerf.o	:$(SCMT)/MixtPerf.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/intpdiff.o	:$(SCMT)/intpdiff.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/inflow_bc.o	:$(SCMT)/inflow_bc.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/outflow_bc.o	:$(SCMT)/outflow_bc.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/wall_bc.o	:$(SCMT)/wall_bc.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/filters_cmt.o	:$(SCMT)/filters_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/diagnostics.o	:$(SCMT)/diagnostics.f;		$(FC) -c $(FL2) $< -o $@

# MXM       ############################################################################
$(OBJDIR)/mxm_wrapper.o	  :$S/mxm_wrapper.f;		$(FC) -c $(FL2) $< -o $@ 
$(OBJDIR)/mxm_std.o	  :$S/mxm_std.f;		$(FC) -c $(FL3) $< -o $@
$(OBJDIR)/mxm_bgq.o	  :$S/mxm_bgq.f;		$(FC) -c $(FL3) $< -o $@

# C Files ##################################################################################
$(OBJDIR)/nek_comm.o             :$S/nek_comm.c;          $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/byte.o                 :$S/byte.c;              $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/chelpers.o             :$S/chelpers.c;          $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)fail.o            :$(J)/fail.c;            $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)tensor.o          :$(J)/tensor.c;          $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)sort.o            :$(J)/sort.c;            $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)sarray_sort.o     :$(J)/sarray_sort.c;     $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)comm.o            :$(J)/comm.c;            $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)crystal.o         :$(J)/crystal.c;         $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)sarray_transfer.o :$(J)/sarray_transfer.c; $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)fcrystal.o        :$(J)/fcrystal.c;        $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)gs.o              :$(J)/gs.c;              $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)gs_local.o        :$(J)/gs_local.c;        $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)fcrs.o            :$(J)/fcrs.c;            $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)xxt.o             :$(J)/xxt.c;             $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)sparse_cholesky.o :$(J)/sparse_cholesky.c; $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)amg.o             :$(J)/amg.c;             $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)findpts.o         :$(J)/findpts.c;         $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)findpts_local.o   :$(J)/findpts_local.c;   $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)obbox.o           :$(J)/obbox.c;           $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)poly.o            :$(J)/poly.c;            $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)lob_bnd.o         :$(J)/lob_bnd.c;         $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)findpts_el_3.o    :$(J)/findpts_el_3.c;    $(CC) -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/$(JO)findpts_el_2.o    :$(J)/findpts_el_2.c;    $(CC) -c $(cFL2) $(JL) $< -o $@

# 3rd party #######################################################################################
$(OBJDIR)/dsygv.o     	:$S/3rd_party/dsygv.f;  	$(FC) -c $(L0) $(PPS_F) $< -o $@
$(OBJDIR)/ssygv.o     	:$S/3rd_party/ssygv.f;  	$(FC) -c $(L0) $(PPS_F) $< -o $@
$(OBJDIR)/blas.o   	:$S/3rd_party/blas.f;		$(FC) -c $(L0) $(PPS_F) $< -o $@
$(OBJDIR)/visit.o        :$S/3rd_party/visit.f;		$(FC) -c $(VISIT_IFLAGS) $(FL2) $< -o $@
$(OBJDIR)/nek_in_situ.o  :$S/3rd_party/nek_in_situ.f;	$(FC) -c $(VISIT_IFLAGS) $(FL2) $< -o $@
$(OBJDIR)/finiparser.o   :$S/3rd_party/finiparser.c;    $(CC)  -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/iniparser.o    :$S/3rd_party/iniparser.c;     $(CC)  -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/dictionary.o   :$S/3rd_party/dictionary.c;    $(CC)  -c $(cFL2) $(JL) $< -o $@
