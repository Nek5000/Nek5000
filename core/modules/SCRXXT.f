      module scrxxt_mod

      implicit none

      real, allocatable, target :: cb_scrxxt(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        navier8.f (set_up_h1_crs/crs_solve_h1: cmlt,mask(lcr,lelv))
         integer sz
         integer lxc, lcr

         lxc = 2
         lcr = lxc**ldim

         sz = 2*lcr*lelv

         n = sz

         allocate(cb_scrxxt(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrxxt$',ierr)
         cb_scrxxt = 0

      end subroutine init
      end module scrxxt_mod
