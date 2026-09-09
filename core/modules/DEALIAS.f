      module dealias_mod

      implicit none
c
c     Dealiasing variables
c
      real, allocatable, target :: cb_solnd(:)
      real, allocatable, target :: cb_interpd(:)

      real, pointer :: vxd(:,:,:,:), vyd(:,:,:,:), vzd(:,:,:,:)

      real, pointer :: imd1(:,:), imd1t(:,:)
      real, pointer :: im1d(:,:), im1dt(:,:)
      real, pointer :: pmd1(:,:), pmd1t(:,:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_solnd(lxd*lyd*lzd*lelv * 3), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_solnd$',ierr)
         cb_solnd = 0

         allocate(cb_interpd(lx1*lxd * 6), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_interpd$',ierr)
         cb_interpd = 0

c        Group 1: /solnd/
         ioff = 1
         vxd(1:lxd,1:lyd,1:lzd,1:lelv) =>
     $         cb_solnd(ioff : ioff + lxd*lyd*lzd*lelv - 1)
         ioff = ioff + lxd*lyd*lzd*lelv
         vyd(1:lxd,1:lyd,1:lzd,1:lelv) =>
     $         cb_solnd(ioff : ioff + lxd*lyd*lzd*lelv - 1)
         ioff = ioff + lxd*lyd*lzd*lelv
         vzd(1:lxd,1:lyd,1:lzd,1:lelv) =>
     $         cb_solnd(ioff : ioff + lxd*lyd*lzd*lelv - 1)

c        Group 2: /interpd/
         ioff = 1
         imd1(1:lx1,1:lxd) =>
     $         cb_interpd(ioff : ioff + lx1*lxd - 1)
         ioff = ioff + lx1*lxd
         imd1t(1:lxd,1:lx1) =>
     $         cb_interpd(ioff : ioff + lx1*lxd - 1)
         ioff = ioff + lx1*lxd
         im1d(1:lxd,1:lx1) =>
     $         cb_interpd(ioff : ioff + lx1*lxd - 1)
         ioff = ioff + lx1*lxd
         im1dt(1:lx1,1:lxd) =>
     $         cb_interpd(ioff : ioff + lx1*lxd - 1)
         ioff = ioff + lx1*lxd
         pmd1(1:lx1,1:lxd) =>
     $         cb_interpd(ioff : ioff + lx1*lxd - 1)
         ioff = ioff + lx1*lxd
         pmd1t(1:lxd,1:lx1) =>
     $         cb_interpd(ioff : ioff + lx1*lxd - 1)

      end subroutine init
      end module dealias_mod
