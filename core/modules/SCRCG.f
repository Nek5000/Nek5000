      module scrcg_mod

      implicit none

      real, allocatable, target :: cb_scrcg(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrcg

c        hmholtz.f (hmh_flex_cg: r,z,p,w, each LX1*LY1*LZ1*LELT)
c        comm_mpi.f gop_test uses nwd(500)
         nscrcg = max(4*lx1*ly1*lz1*lelt, 500)

         allocate(cb_scrcg(nscrcg), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrcg$',ierr)
         cb_scrcg = 0

      end subroutine init
      end module scrcg_mod
