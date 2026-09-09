      module vproj_mod

      implicit none
c
c     Vproj variables
c
      integer, allocatable, target :: cb_wrthoi(:)
      real,    allocatable, target :: cb_wrthov(:)

      integer, pointer :: ivproj(:,:)

      real, pointer :: vproj(:,:)   ! independent storage
                                     ! (formerly overlapped ORTHOSTRS)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff, ktop

         ktop = lelv*lx1*ly1*lz1

c        --- allocate backing arrays ---

         allocate(cb_wrthoi(2*6), stat=ierr)             ! ivproj
         if (ierr.ne.0) call exitti('alloc cb_wrthoi$',ierr)
         cb_wrthoi = 0

         allocate(cb_wrthov(ktop*mxprev*2*ldim), stat=ierr) ! vproj
         if (ierr.ne.0) call exitti('alloc cb_wrthov$',ierr)
         cb_wrthov = 0

c        Group 1: /wrthoi/
         ioff = 1
         ivproj(1:2,1:6) => cb_wrthoi(ioff : ioff + 2*6 - 1)

c        Group 2: /wrthov/
         ioff = 1
         vproj(1:ktop*mxprev,1:2*ldim) =>
     $         cb_wrthov(ioff : ioff + ktop*mxprev*2*ldim - 1)

      end subroutine init
      end module vproj_mod
