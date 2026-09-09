      module scrch_mod

      implicit none

      real, allocatable, target :: cb_scrch(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrch

c        navier6.f (z(2*ltotd)); dominates every other site since
c        those are all single LX1*LY1*LZ1*LELT/LELV-sized arrays and
c        LELV<=LELT
         nscrch = 2*lx1*ly1*lz1*lelt

         allocate(cb_scrch(nscrch), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrch$',ierr)
         cb_scrch = 0

      end subroutine init
      end module scrch_mod
