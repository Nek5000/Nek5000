AUTOMAKE_OPTIONS = foreign subdir-objects

include ../Makefile.problem

PROBLEM         = turbChannel
DRIVER          = single_physics_turbChannel

single_physics_turbChannel_SOURCES = ../../../../tests/drivers/single_physics.cxx
single_physics_turbChannel_LDADD = $(LDADD)
single_physics_turbChannel_CPPFLAGS = -DENABLE_NEK

