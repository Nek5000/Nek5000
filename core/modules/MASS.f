      module mass_mod

      implicit none
c
c     Mass matrix
c
      real, allocatable, target :: cb_mass(:)

      real, pointer :: bm1(:,:,:,:), bm2(:,:,:,:)
      real, pointer :: binvm1(:,:,:,:), bintm1(:,:,:,:)
      real, pointer :: bm2inv(:,:,:,:), baxm1(:,:,:,:)
      real, pointer :: bm1lag(:,:,:,:,:)
      real, pointer :: volvm1, volvm2, voltm1, voltm2
      real, pointer :: yinvm1(:,:,:,:)
      real, pointer :: binvdg(:,:)
      real, pointer :: bm1ms(:,:,:,:)   !weighted mass matrix
      real, pointer :: upf(:,:,:,:)     !unity partition function
      real, pointer :: volvm1ms

      contains

      subroutine init
         use size_mod
         implicit none
         integer ltot, ltotv, ltot2

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_mass(
     $        lx1*ly1*lz1*lelt              ! bm1
     $      + lx2*ly2*lz2*lelv              ! bm2
     $      + lx1*ly1*lz1*lelv              ! binvm1
     $      + lx1*ly1*lz1*lelt              ! bintm1
     $      + lx2*ly2*lz2*lelt              ! bm2inv
     $      + lx1*ly1*lz1*lelt              ! baxm1
     $      + lx1*ly1*lz1*lelt*(lorder-1)   ! bm1lag
     $      + 1                             ! volvm1
     $      + 1                             ! volvm2
     $      + 1                             ! voltm1
     $      + 1                             ! voltm2
     $      + lx1*ly1*lz1*lelt              ! yinvm1
     $      + lx1*ly1*lz1*lelt              ! binvdg
     $      + lx1*ly1*lz1*lelt              ! bm1ms
     $      + lx1*ly1*lz1*lelt              ! upf
     $      + 1), stat=ierr)                ! volvm1ms
         if (ierr.ne.0) call exitti('alloc cb_mass$',ierr)
         cb_mass = 0

c        Group 1: /mass/
         ioff = 1
         ltot = lx1*ly1*lz1*lelt
         bm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         ltot2 = lx2*ly2*lz2*lelv
         bm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_mass(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         ltotv = lx1*ly1*lz1*lelv
         binvm1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_mass(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         bintm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         bm2inv(1:lx2,1:ly2,1:lz2,1:lelt) =>
     $         cb_mass(ioff : ioff + lx2*ly2*lz2*lelt - 1)
         ioff = ioff + lx2*ly2*lz2*lelt
         baxm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         bm1lag(1:lx1,1:ly1,1:lz1,1:lelt,1:lorder-1) =>
     $         cb_mass(ioff : ioff + lx1*ly1*lz1*lelt*(lorder-1) - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*(lorder-1)
         volvm1 => cb_mass(ioff)
         ioff = ioff + 1
         volvm2 => cb_mass(ioff)
         ioff = ioff + 1
         voltm1 => cb_mass(ioff)
         ioff = ioff + 1
         voltm2 => cb_mass(ioff)
         ioff = ioff + 1
         yinvm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         binvdg(1:lx1*ly1*lz1,1:lelt) => cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         bm1ms(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         upf(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_mass(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         volvm1ms => cb_mass(ioff)

      end subroutine init
      end module mass_mod
