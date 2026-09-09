      module scrpre_mod

      implicit none

      real, allocatable, target :: cb_scrpre(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        fasts.f (local_solves_fdm: v1(lx1*ly1*lz1*lelv),
c        w1,w2(lx1*ly1*lz1)); dominates navier8.f's crs_solve_l2
c        (uc(lcr*lelt),w(2*lx1*ly1*lz1)) and crs_solve_h1 (uc(lcr*lelt))
         integer sz

         sz = lx1*ly1*lz1*lelv + 2*lx1*ly1*lz1

         n = sz

         allocate(cb_scrpre(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrpre$',ierr)
         cb_scrpre = 0

      end subroutine init
      end module scrpre_mod
