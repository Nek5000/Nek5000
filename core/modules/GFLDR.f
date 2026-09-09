      module gfldr_mod

      include 'SIZE'
      implicit none

      real             tol, bb_t
      parameter(tol=5e-13, bb_t=0.01)
      integer          ltots
      parameter(ltots=4*lx1*ly1*lz1*lelt) ! max source mesh size

c     lazy allocation -- gfldr() is an optional, rarely-used feature,
c     so the (potentially large) buffers below are only allocated the
c     first time gfldr() is actually called, not at program startup
      logical isinit
      save    isinit
      data    isinit /.false./

      integer                fldh_gfldr, inth_gfldr
      logical                ifbswp
      integer*8              noff0_b, ntots_b, rankoff_b, nSizeFld_b
      integer                nelgs, nels, nxyzs

      real,    allocatable, target :: cb_gfldr_r(:)
      real*4,  allocatable, target :: cb_gfldr_r4(:)
      integer, allocatable, target :: cb_gfldr_i(:)

      real,    pointer :: buffld(:)  ! one source field
      real,    pointer :: xm1s(:)    ! source x coord
      real,    pointer :: ym1s(:)    ! source y coord
      real,    pointer :: zm1s(:)    ! source z coord
      real,    pointer :: grst(:)    ! interpolation r-s-t working array
      real,    pointer :: gdist(:)

      real*4,  pointer :: bufr(:)    ! read buffer

      integer, pointer :: grcode(:)
      integer, pointer :: gelid(:)
      integer, pointer :: gproc(:)

      contains

      subroutine init
         implicit none

         integer ierr, ioff

         if (isinit) return
         isinit = .true.

c        --- allocate backing arrays ---

         allocate(cb_gfldr_r(ltots*(5+ldim)), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_gfldr_r$',ierr)
         cb_gfldr_r = 0

         allocate(cb_gfldr_r4(2*ldim*ltots), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_gfldr_r4$',ierr)
         cb_gfldr_r4 = 0

         allocate(cb_gfldr_i(3*ltots), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_gfldr_i$',ierr)
         cb_gfldr_i = 0

c        Group 1: real work arrays
         ioff = 1
         buffld(1:ltots) => cb_gfldr_r(ioff : ioff + ltots - 1)
         ioff = ioff + ltots
         xm1s(1:ltots) => cb_gfldr_r(ioff : ioff + ltots - 1)
         ioff = ioff + ltots
         ym1s(1:ltots) => cb_gfldr_r(ioff : ioff + ltots - 1)
         ioff = ioff + ltots
         zm1s(1:ltots) => cb_gfldr_r(ioff : ioff + ltots - 1)
         ioff = ioff + ltots
         grst(1:ldim*ltots) => cb_gfldr_r(ioff : ioff + ldim*ltots - 1)
         ioff = ioff + ldim*ltots
         gdist(1:ltots) => cb_gfldr_r(ioff : ioff + ltots - 1)

c        Group 2: real*4 read buffer
         bufr(1:2*ldim*ltots) => cb_gfldr_r4(1 : 2*ldim*ltots)

c        Group 3: integer work arrays
         ioff = 1
         grcode(1:ltots) => cb_gfldr_i(ioff : ioff + ltots - 1)
         ioff = ioff + ltots
         gelid(1:ltots) => cb_gfldr_i(ioff : ioff + ltots - 1)
         ioff = ioff + ltots
         gproc(1:ltots) => cb_gfldr_i(ioff : ioff + ltots - 1)

      end subroutine init
      end module gfldr_mod
