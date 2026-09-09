      module orthov_mod

      implicit none

      real, allocatable, target :: cb_orthov(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, northov
         integer nset

c        induct.f (pset: LX2*LY2*LZ2*LELV*MXPREV*NSET); this also
c        covers navier4.f's RHS(LTOT2,MXPREV) site since NSET>=1
         nset = 1 + lbelv/lelv

         northov = lx2*ly2*lz2*lelv*mxprev*nset

         allocate(cb_orthov(northov), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_orthov$',ierr)
         cb_orthov = 0

      end subroutine init
      end module orthov_mod
