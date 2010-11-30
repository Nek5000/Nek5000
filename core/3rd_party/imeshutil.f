c random utilities to help with iMesh
c called by MYASSERT macro
#define NULLSTRIP(s) s(:index(s, char(0))-1)

c building without imesh for now

      subroutine imesh_err(ierr, imesh, file, line)

      print *, "ASSERT ERROR: ", ierr, 'in ', file, ' line', line

#if 0
  implicit none
  integer ierr, line
  integer*8 imesh
  character file(*)
  character*1024 errmsg



  call iMesh_getDescription(imesh, &!iMesh_Instance instance,
                            errmsg, & !/*inout*/ char *descr,
                            ierr)
  print *, NULLSTRIP(errmsg)

!  call exit
#endif
      return
      end 


