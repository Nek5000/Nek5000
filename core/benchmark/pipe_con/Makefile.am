AUTOMAKE_OPTIONS = foreign subdir-objects

include ../Makefile.problem

PROBLEM         = pipe_con
DRIVER          = single_physics_pipe_con

single_physics_pipe_con_SOURCES = ../../../../tests/drivers/single_physics.cxx
single_physics_pipe_con_LDADD = $(LDADD)
single_physics_pipe_con_CPPFLAGS = -DENABLE_NEK

