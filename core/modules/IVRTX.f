      module ivrtx_mod

      implicit none

      integer*8, allocatable, target :: cb_ivrtx(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        vertex((2**ldim)*lelt), integer*8; shared across connect1.f,
c        conduct.f, drive1.f, gmres.f, hrefine.f, hsmg.f, map2.f,
c        navier0.f, navier8.f, subs1.f (13 sites)
         integer sz

         sz = (2**ldim)*lelt

         n = sz

         allocate(cb_ivrtx(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_ivrtx$',ierr)
         cb_ivrtx = 0_8

      end subroutine init
      end module ivrtx_mod
