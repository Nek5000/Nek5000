      module cvflow_a_mod

      implicit none

      real, allocatable, target :: cb_cvflow_a(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        drive2.f/multimesh.f (vol_flow/vol_flow_ms:
c        vxc,vyc,vzc(lx1*ly1*lz1*lelv), prc(lx2*ly2*lz2*lelv),
c        vdc(lx1*ly1*lz1*lelv,2))
         integer sz

         sz = 5*lx1*ly1*lz1*lelv + lx2*ly2*lz2*lelv

         n = sz

         allocate(cb_cvflow_a(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_cvflow_a$',ierr)
         cb_cvflow_a = 0

      end subroutine init
      end module cvflow_a_mod
