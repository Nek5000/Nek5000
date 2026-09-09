      module fdmh1_mod
      include 'SIZE'
c
c     'FDMH1'
c
      real, allocatable, target :: cb_fdmh1r(:)

      real, pointer :: dd(:,:), fds(:,:), fdst(:,:)
      real, pointer :: elsize(:,:), bhalf(:,:,:,:)
c
      integer         kfldfdm
      common /fdmh1i/ kfldfdm

      integer, allocatable, target :: ktype(:,:,:)
c
      logical         ifbhalf
      common /fdmh1l/ ifbhalf

      contains

      subroutine init
         implicit none
         integer ierr, ioff

         allocate(cb_fdmh1r(9*lx1 + 2*9*lx1*lx1 + 3*lelt
     $                     + lx1*ly1*lz1*lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_fdmh1r$',ierr)
         cb_fdmh1r = 0
         allocate(ktype(lelt,3,0:4), stat=ierr)
         if (ierr.ne.0) call exitti('alloc ktype$',ierr)
         ktype = 0

         ioff = 1
         dd(1:lx1,1:9) => cb_fdmh1r(ioff : ioff+lx1*9-1)
         ioff = ioff + lx1*9
         fds(1:lx1*lx1,1:9) => cb_fdmh1r(ioff : ioff+lx1*lx1*9-1)
         ioff = ioff + lx1*lx1*9
         fdst(1:lx1*lx1,1:9) => cb_fdmh1r(ioff : ioff+lx1*lx1*9-1)
         ioff = ioff + lx1*lx1*9
         elsize(1:3,1:lelt) => cb_fdmh1r(ioff : ioff+3*lelt-1)
         ioff = ioff + 3*lelt
         bhalf(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $      cb_fdmh1r(ioff : ioff+lx1*ly1*lz1*lelt-1)

      end subroutine init
      end module fdmh1_mod
