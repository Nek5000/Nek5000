      module mvgeom_mod

      implicit none
c
c     Moving mesh variables
c
      real, allocatable, target :: cb_wsol(:)
      real, allocatable, target :: cb_wlag(:)
      real, allocatable, target :: cb_wmsu(:)
      real, allocatable, target :: cb_eigvec(:)

      real, pointer :: wx(:,:,:,:), wy(:,:,:,:), wz(:,:,:,:)

      real, pointer :: wxlag(:,:,:,:,:), wylag(:,:,:,:,:)
     $               , wzlag(:,:,:,:,:)

      real, pointer :: w1mask(:,:,:,:), w2mask(:,:,:,:)
     $               , w3mask(:,:,:,:), wmult(:,:,:,:)

      real, pointer :: ev1(:,:,:,:), ev2(:,:,:,:), ev3(:,:,:,:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

c        --- allocate backing arrays ---

         allocate(cb_wsol(lx1m*ly1m*lz1m*lelt * 3), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_wsol$',ierr)
         cb_wsol = 0

         allocate(cb_wlag(lx1m*ly1m*lz1m*lelt*(lorder-1) * 3),
     $            stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_wlag$',ierr)
         cb_wlag = 0

         allocate(cb_wmsu(lx1m*ly1m*lz1m*lelt * 4), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_wmsu$',ierr)
         cb_wmsu = 0

         allocate(cb_eigvec(lx1m*ly1m*lz1m*lelv * 3), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_eigvec$',ierr)
         cb_eigvec = 0

c        Group 1: /wsol/
         ioff = 1
         wx(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wsol(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         wy(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wsol(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         wz(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wsol(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)

c        Group 2: /wlag/
         ioff = 1
         wxlag(1:lx1m,1:ly1m,1:lz1m,1:lelt,1:lorder-1) =>
     $         cb_wlag(ioff : ioff
     $         + lx1m*ly1m*lz1m*lelt*(lorder-1) - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt*(lorder-1)
         wylag(1:lx1m,1:ly1m,1:lz1m,1:lelt,1:lorder-1) =>
     $         cb_wlag(ioff : ioff
     $         + lx1m*ly1m*lz1m*lelt*(lorder-1) - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt*(lorder-1)
         wzlag(1:lx1m,1:ly1m,1:lz1m,1:lelt,1:lorder-1) =>
     $         cb_wlag(ioff : ioff
     $         + lx1m*ly1m*lz1m*lelt*(lorder-1) - 1)

c        Group 3: /wmsu/
         ioff = 1
         w1mask(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wmsu(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         w2mask(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wmsu(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         w3mask(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wmsu(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         wmult(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_wmsu(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)

c        Group 4: /eigvec/
         ioff = 1
         ev1(1:lx1m,1:ly1m,1:lz1m,1:lelv) =>
     $         cb_eigvec(ioff : ioff + lx1m*ly1m*lz1m*lelv - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelv
         ev2(1:lx1m,1:ly1m,1:lz1m,1:lelv) =>
     $         cb_eigvec(ioff : ioff + lx1m*ly1m*lz1m*lelv - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelv
         ev3(1:lx1m,1:ly1m,1:lz1m,1:lelv) =>
     $         cb_eigvec(ioff : ioff + lx1m*ly1m*lz1m*lelv - 1)

      end subroutine init
      end module mvgeom_mod
