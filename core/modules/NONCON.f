      module noncon_mod

      implicit none
c
c     Non conforming variables
c
      real,    allocatable, target :: cb_allr(:)
      integer, allocatable, target :: cb_alli(:)
      logical, allocatable, target :: cb_logg(:)

      real, pointer :: umult(:)
      real, pointer :: jmat(:,:,:,:)
      real, pointer :: xsp(:,:), ysp(:,:), zsp(:,:)
      real, pointer :: xch(:,:), ych(:,:), zch(:,:)
      real, pointer :: rtwid(:), stwid(:)
      real, pointer :: dtrk(:)
      real, pointer :: rs(:,:,:)

      integer, pointer :: noncon_f(:)
      integer, pointer :: noncon_e(:)
      integer, pointer :: noncon_ip(:)
      integer, pointer :: mortar(:,:)
      integer, pointer :: imin(:,:)
      integer, pointer :: mort_m

      logical, pointer :: ifnc, ifhalf
      logical, pointer :: ifjt(:)

      contains

      subroutine init
         use size_mod
         implicit none
         integer ltot

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_allr(
     $        lx1*ly1*lz1*lelt          ! umult
     $      + lx1*lx1*2*maxmor          ! jmat
     $      + lx1*lz1                   ! xsp
     $      + lx1*lz1                   ! ysp
     $      + lx1*lz1                   ! zsp
     $      + lx1*lz1                   ! xch
     $      + lx1*lz1                   ! ych
     $      + lx1*lz1                   ! zch
     $      + lx1                       ! rtwid
     $      + lx1                       ! stwid
     $      + ldim                      ! dtrk
     $      + 2*2*2), stat=ierr)        ! rs
         if (ierr.ne.0) call exitti('alloc cb_allr$',ierr)
         cb_allr = 0

         allocate(cb_alli(
     $        maxmor                    ! noncon_f
     $      + maxmor                    ! noncon_e
     $      + maxmor                    ! noncon_ip
     $      + 6*lelt                    ! mortar
     $      + 3*2                       ! imin
     $      + 1), stat=ierr)            ! mort_m
         if (ierr.ne.0) call exitti('alloc cb_alli$',ierr)
         cb_alli = 0

         allocate(cb_logg(
     $        1                         ! ifnc
     $      + 1                         ! ifhalf
     $      + maxmor), stat=ierr)       ! ifjt
         if (ierr.ne.0) call exitti('alloc cb_logg$',ierr)
         cb_logg = .false.

c        Group 1: /allr/
         ioff = 1
         ltot = lx1*ly1*lz1*lelt
         umult(1:lx1*ly1*lz1*lelt) => cb_allr(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         jmat(1:lx1,1:lx1,1:2,1:maxmor) =>
     $         cb_allr(ioff : ioff + lx1*lx1*2*maxmor - 1)
         ioff = ioff + lx1*lx1*2*maxmor
         xsp(1:lx1,1:lz1) => cb_allr(ioff : ioff + lx1*lz1 - 1)
         ioff = ioff + lx1*lz1
         ysp(1:lx1,1:lz1) => cb_allr(ioff : ioff + lx1*lz1 - 1)
         ioff = ioff + lx1*lz1
         zsp(1:lx1,1:lz1) => cb_allr(ioff : ioff + lx1*lz1 - 1)
         ioff = ioff + lx1*lz1
         xch(1:lx1,1:lz1) => cb_allr(ioff : ioff + lx1*lz1 - 1)
         ioff = ioff + lx1*lz1
         ych(1:lx1,1:lz1) => cb_allr(ioff : ioff + lx1*lz1 - 1)
         ioff = ioff + lx1*lz1
         zch(1:lx1,1:lz1) => cb_allr(ioff : ioff + lx1*lz1 - 1)
         ioff = ioff + lx1*lz1
         rtwid(1:lx1) => cb_allr(ioff : ioff + lx1 - 1)
         ioff = ioff + lx1
         stwid(1:lx1) => cb_allr(ioff : ioff + lx1 - 1)
         ioff = ioff + lx1
         dtrk(1:ldim) => cb_allr(ioff : ioff + ldim - 1)
         ioff = ioff + ldim
         rs(1:2,1:2,1:2) => cb_allr(ioff : ioff + 2*2*2 - 1)

c        Group 2: /alli/
         ioff = 1
         noncon_f(1:maxmor) => cb_alli(ioff : ioff + maxmor - 1)
         ioff = ioff + maxmor
         noncon_e(1:maxmor) => cb_alli(ioff : ioff + maxmor - 1)
         ioff = ioff + maxmor
         noncon_ip(1:maxmor) => cb_alli(ioff : ioff + maxmor - 1)
         ioff = ioff + maxmor
         mortar(1:6,1:lelt) => cb_alli(ioff : ioff + 6*lelt - 1)
         ioff = ioff + 6*lelt
         imin(1:3,1:2) => cb_alli(ioff : ioff + 3*2 - 1)
         ioff = ioff + 3*2
         mort_m => cb_alli(ioff)

c        Group 3: /logg/
         ioff = 1
         ifnc => cb_logg(ioff)
         ioff = ioff + 1
         ifhalf => cb_logg(ioff)
         ioff = ioff + 1
         ifjt(1:maxmor) => cb_logg(ioff : ioff + maxmor - 1)

      end subroutine init
      end module noncon_mod
