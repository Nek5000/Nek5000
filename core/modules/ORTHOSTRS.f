      module orthostrs_mod

      implicit none

      real,    allocatable, target :: cb_wrthov(:)
      integer, allocatable, target :: cb_srthoi(:)

      real, pointer :: xstrs(:)   ! independent storage
      real, pointer :: bstrs(:)   ! (formerly overlapped VPROJ)

      integer, pointer :: napproxstrs(:)

      character(len=4) :: name4strs

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff, ktotstrs

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

         ktotstrs = lx1m*ly1m*lz1m*lelt

c        --- allocate backing arrays ---

         allocate(cb_wrthov(
     $        ldim*ktotstrs*(1+mxprev)     ! xstrs
     $      + ldim*ktotstrs*(1+mxprev)),   ! bstrs
     $      stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_wrthov$',ierr)
         cb_wrthov = 0

         allocate(cb_srthoi(2), stat=ierr)  ! napproxstrs
         if (ierr.ne.0) call exitti('alloc cb_srthoi$',ierr)
         cb_srthoi = 0

c        Group 1: /wrthov/
         ioff = 1
         xstrs(1:ldim*ktotstrs*(1+mxprev)) =>
     $         cb_wrthov(ioff : ioff + ldim*ktotstrs*(1+mxprev) - 1)
         ioff = ioff + ldim*ktotstrs*(1+mxprev)
         bstrs(1:ldim*ktotstrs*(1+mxprev)) =>
     $         cb_wrthov(ioff : ioff + ldim*ktotstrs*(1+mxprev) - 1)

c        Group 2: /srthoi/
         ioff = 1
         napproxstrs(1:2) => cb_srthoi(ioff : ioff + 2 - 1)

      end subroutine init
      end module orthostrs_mod
