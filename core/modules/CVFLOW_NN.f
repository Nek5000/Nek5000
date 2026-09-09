      module cvflow_nn_mod

      implicit none

      real, allocatable, target :: cb_cvflow_nn(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        multimesh.f (plan3_vol_ms/plan4_vol_ms:
c        vxcbc,vycbc,vzcbc(lx1*ly1*lz1*lelv))
         integer sz

         sz = 3*lx1*ly1*lz1*lelv

         n = sz

         allocate(cb_cvflow_nn(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_cvflow_nn$',ierr)
         cb_cvflow_nn = 0

      end subroutine init
      end module cvflow_nn_mod
