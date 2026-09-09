      module scrpr2_mod

      implicit none

      real, allocatable, target :: cb_scrpr2(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        subs1.f (crs_strs: vc1,vc2,vc3(lcr*lelt)); dominates
c        navier8.f (crs_solve_h1: vc(lcr*lelt))
         integer sz
         integer lxc, lcr

         lxc = 2
         lcr = lxc**ldim

         sz = 3*lcr*lelt

         n = sz

         allocate(cb_scrpr2(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrpr2$',ierr)
         cb_scrpr2 = 0

      end subroutine init
      end module scrpr2_mod
