      module swaplengths_mod

      implicit none

      real, allocatable, target :: cb_swaplengths(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        fast3d.f (swap_lengths) / fasts.f (init_weight_op):
c        l(lx1,ly1,lz1,lelv)
         integer sz

         sz = lx1*ly1*lz1*lelv

         n = sz

         allocate(cb_swaplengths(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_swaplengths$',ierr)
         cb_swaplengths = 0

      end subroutine init
      end module swaplengths_mod
