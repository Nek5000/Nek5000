      module weightop_mod

      implicit none

      real, allocatable, target :: cb_weightop(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        fasts.f (init_weight_op/do_weight_op): w(lx2,lz2,2,3,levb)
         integer sz
         integer levb

         levb = lelv + lbelv

         sz = lx2*lz2*2*3*levb

         n = sz

         allocate(cb_weightop(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_weightop$',ierr)
         cb_weightop = 0

      end subroutine init
      end module weightop_mod
