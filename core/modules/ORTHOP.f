      module orthop_mod

      implicit none

      integer laxtp   ! max dim of projection space (depends on SIZE)

      integer, allocatable, target :: cb_prthoi(:)
      real,    allocatable, target :: cb_prthov(:)

      integer, pointer :: napproxp(:)

      real, pointer :: approxp(:,:)

      character(len=4) :: name4p

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff, ktotp

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

         laxtp = mxprev
         ktotp = lx1*ly1*lz1*lelt

c        --- allocate backing arrays ---

         allocate(cb_prthoi(2), stat=ierr)             ! napproxp
         if (ierr.ne.0) call exitti('alloc cb_prthoi$',ierr)
         cb_prthoi = 0

         allocate(cb_prthov(ktotp*(laxtp+1)), stat=ierr) ! approxp
         if (ierr.ne.0) call exitti('alloc cb_prthov$',ierr)
         cb_prthov = 0

c        Group 1: /prthoi/
         ioff = 1
         napproxp(1:2) => cb_prthoi(ioff : ioff + 2 - 1)

c        Group 2: /prthov/
         ioff = 1
         approxp(1:ktotp,0:laxtp) =>
     $         cb_prthov(ioff : ioff + ktotp*(laxtp+1) - 1)

      end subroutine init
      end module orthop_mod
