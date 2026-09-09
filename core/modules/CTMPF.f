      module ctmpf_mod

      implicit none

      real, allocatable, target :: cb_ctmpf(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        fast3d.f (gen_fast/gen_fast_spacing/swap_lengths):
c        lr,ls,lt(2*lx1+4), llr,lls,llt,lmr,lms,lmt,lrr,lrs,lrt(lelt)
         integer sz

         sz = 3*(2*lx1+4) + 9*lelt

         n = sz

         allocate(cb_ctmpf(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_ctmpf$',ierr)
         cb_ctmpf = 0

      end subroutine init
      end module ctmpf_mod
