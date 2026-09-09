      module orthot_mod

      implicit none

      integer laxtt   ! max dim of projection space (depends on SIZE)

      integer, allocatable, target :: cb_trthoi(:)
      real,    allocatable, target :: cb_trthov(:)

      integer, pointer :: napproxt(:,:)

      real, pointer :: approxt(:,:,:)

      character(len=4) :: name4t

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff, ktott

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

         laxtt = mxprev
         ktott = lx1*ly1*lz1*lelt

c        --- allocate backing arrays ---

         allocate(cb_trthoi(2*ldimt_proj), stat=ierr) ! napproxt
         if (ierr.ne.0) call exitti('alloc cb_trthoi$',ierr)
         cb_trthoi = 0

         allocate(cb_trthov(ktott*(laxtt+1)*ldimt_proj), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_trthov$',ierr)
         cb_trthov = 0
                                                        ! approxt

c        Group 1: /trthoi/
         ioff = 1
         napproxt(1:2,1:ldimt_proj) =>
     $         cb_trthoi(ioff : ioff + 2*ldimt_proj - 1)

c        Group 2: /trthov/
         ioff = 1
         approxt(1:ktott,0:laxtt,1:ldimt_proj) =>
     $         cb_trthov(ioff : ioff + ktott*(laxtt+1)*ldimt_proj - 1)

      end subroutine init
      end module orthot_mod
