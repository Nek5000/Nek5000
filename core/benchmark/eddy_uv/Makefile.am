AUTOMAKE_OPTIONS = foreign subdir-objects

include ../Makefile.problem

PROBLEM         = eddy_uv
DRIVER          = single_physics_eddy_uv

single_physics_eddy_uv_SOURCES = ../../../../tests/drivers/single_physics.cxx
single_physics_eddy_uv_LDADD = $(LDADD)
single_physics_eddy_uv_CPPFLAGS = -DENABLE_NEK

