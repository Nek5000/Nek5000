      module geom_mod

      implicit none
c
c     Geometry arrays
c
      real,    allocatable, target :: cb_gxyz(:)
      real,    allocatable, target :: cb_giso1(:)
      real,    allocatable, target :: cb_giso2(:)
      real,    allocatable, target :: cb_gisod(:)
      real,    allocatable, target :: cb_gmfact(:)
      real,    allocatable, target :: cb_gsurf(:)
      real,    allocatable, target :: cb_gvolm(:)
      logical, allocatable, target :: cb_glog(:)
      integer, allocatable, target :: cb_cbbid(:)

      real, pointer :: xm1(:,:,:,:), ym1(:,:,:,:), zm1(:,:,:,:)
      real, pointer :: xm2(:,:,:,:), ym2(:,:,:,:), zm2(:,:,:,:)

      real, pointer :: rxm1(:,:,:,:), sxm1(:,:,:,:), txm1(:,:,:,:)
      real, pointer :: rym1(:,:,:,:), sym1(:,:,:,:), tym1(:,:,:,:)
      real, pointer :: rzm1(:,:,:,:), szm1(:,:,:,:), tzm1(:,:,:,:)
      real, pointer :: jacm1(:,:,:,:)
      real, pointer :: jacmi(:,:)

      real, pointer :: rxm2(:,:,:,:), sxm2(:,:,:,:), txm2(:,:,:,:)
      real, pointer :: rym2(:,:,:,:), sym2(:,:,:,:), tym2(:,:,:,:)
      real, pointer :: rzm2(:,:,:,:), szm2(:,:,:,:), tzm2(:,:,:,:)
      real, pointer :: jacm2(:,:,:,:)

      real, pointer :: rx(:,:,:)

      real, pointer :: g1m1(:,:,:,:), g2m1(:,:,:,:), g3m1(:,:,:,:)
      real, pointer :: g4m1(:,:,:,:), g5m1(:,:,:,:), g6m1(:,:,:,:)

      real, pointer :: unr(:,:,:), uns(:,:,:), unt(:,:,:)
      real, pointer :: unx(:,:,:,:), uny(:,:,:,:), unz(:,:,:,:)
      real, pointer :: t1x(:,:,:,:), t1y(:,:,:,:), t1z(:,:,:,:)
      real, pointer :: t2x(:,:,:,:), t2y(:,:,:,:), t2z(:,:,:,:)
      real, pointer :: area(:,:,:,:)
      real, pointer :: etalph(:,:,:)
      real, pointer :: dlam

      real, pointer :: vnx(:,:,:,:), vny(:,:,:,:), vnz(:,:,:,:)
      real, pointer :: v1x(:,:,:,:), v1y(:,:,:,:), v1z(:,:,:,:)
      real, pointer :: v2x(:,:,:,:), v2y(:,:,:,:), v2z(:,:,:,:)

      logical, pointer :: ifgeom, ifgmsh3, ifvcor, ifsurt
      logical, pointer :: ifmelt, ifwcno
      logical, pointer :: ifrzer(:)
      logical, pointer :: ifqinp(:,:), ifeppm(:,:)
      logical, pointer :: iflmsf(:), iflmse(:), iflmsc(:)
      logical, pointer :: ifmsfc(:,:,:)
      logical, pointer :: ifmseg(:,:,:)
      logical, pointer :: ifmscr(:,:,:)
      logical, pointer :: ifnskp(:,:)
      logical, pointer :: ifbcor

      integer, pointer :: boundaryID(:,:), boundaryIDt(:,:)

      contains

      subroutine init
         use size_mod
         implicit none
         integer ltot, ltot2

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_gxyz(
     $        lx1*ly1*lz1*lelt         ! xm1
     $      + lx1*ly1*lz1*lelt         ! ym1
     $      + lx1*ly1*lz1*lelt         ! zm1
     $      + lx2*ly2*lz2*lelv         ! xm2
     $      + lx2*ly2*lz2*lelv         ! ym2
     $      + lx2*ly2*lz2*lelv), stat=ierr) ! zm2
         if (ierr.ne.0) call exitti('alloc cb_gxyz$',ierr)
         cb_gxyz = 0

         allocate(cb_giso1(
     $        lx1*ly1*lz1*lelt         ! rxm1
     $      + lx1*ly1*lz1*lelt         ! sxm1
     $      + lx1*ly1*lz1*lelt         ! txm1
     $      + lx1*ly1*lz1*lelt         ! rym1
     $      + lx1*ly1*lz1*lelt         ! sym1
     $      + lx1*ly1*lz1*lelt         ! tym1
     $      + lx1*ly1*lz1*lelt         ! rzm1
     $      + lx1*ly1*lz1*lelt         ! szm1
     $      + lx1*ly1*lz1*lelt         ! tzm1
     $      + lx1*ly1*lz1*lelt         ! jacm1
     $      + lx1*ly1*lz1*lelt), stat=ierr) ! jacmi
         if (ierr.ne.0) call exitti('alloc cb_giso1$',ierr)
         cb_giso1 = 0

         allocate(cb_giso2(
     $        lx2*ly2*lz2*lelv         ! rxm2
     $      + lx2*ly2*lz2*lelv         ! sxm2
     $      + lx2*ly2*lz2*lelv         ! txm2
     $      + lx2*ly2*lz2*lelv         ! rym2
     $      + lx2*ly2*lz2*lelv         ! sym2
     $      + lx2*ly2*lz2*lelv         ! tym2
     $      + lx2*ly2*lz2*lelv         ! rzm2
     $      + lx2*ly2*lz2*lelv         ! szm2
     $      + lx2*ly2*lz2*lelv         ! tzm2
     $      + lx2*ly2*lz2*lelv), stat=ierr) ! jacm2
         if (ierr.ne.0) call exitti('alloc cb_giso2$',ierr)
         cb_giso2 = 0

         allocate(cb_gisod(lxd*lyd*lzd*ldim*ldim*lelv), stat=ierr) ! rx
         if (ierr.ne.0) call exitti('alloc cb_gisod$',ierr)
         cb_gisod = 0

         allocate(cb_gmfact(
     $        lx1*ly1*lz1*lelt         ! g1m1
     $      + lx1*ly1*lz1*lelt         ! g2m1
     $      + lx1*ly1*lz1*lelt         ! g3m1
     $      + lx1*ly1*lz1*lelt         ! g4m1
     $      + lx1*ly1*lz1*lelt         ! g5m1
     $      + lx1*ly1*lz1*lelt), stat=ierr) ! g6m1
         if (ierr.ne.0) call exitti('alloc cb_gmfact$',ierr)
         cb_gmfact = 0

         allocate(cb_gsurf(
     $        lx1*lz1*6*lelt           ! unr
     $      + lx1*lz1*6*lelt           ! uns
     $      + lx1*lz1*6*lelt           ! unt
     $      + lx1*lz1*6*lelt           ! unx
     $      + lx1*lz1*6*lelt           ! uny
     $      + lx1*lz1*6*lelt           ! unz
     $      + lx1*lz1*6*lelt           ! t1x
     $      + lx1*lz1*6*lelt           ! t1y
     $      + lx1*lz1*6*lelt           ! t1z
     $      + lx1*lz1*6*lelt           ! t2x
     $      + lx1*lz1*6*lelt           ! t2y
     $      + lx1*lz1*6*lelt           ! t2z
     $      + lx1*lz1*6*lelt           ! area
     $      + lx1*lz1*2*ldim*lelt      ! etalph
     $      + 1), stat=ierr)           ! dlam
         if (ierr.ne.0) call exitti('alloc cb_gsurf$',ierr)
         cb_gsurf = 0

         allocate(cb_gvolm(
     $        lx1m*ly1m*lz1m*lelt      ! vnx
     $      + lx1m*ly1m*lz1m*lelt      ! vny
     $      + lx1m*ly1m*lz1m*lelt      ! vnz
     $      + lx1m*ly1m*lz1m*lelt      ! v1x
     $      + lx1m*ly1m*lz1m*lelt      ! v1y
     $      + lx1m*ly1m*lz1m*lelt      ! v1z
     $      + lx1m*ly1m*lz1m*lelt      ! v2x
     $      + lx1m*ly1m*lz1m*lelt      ! v2y
     $      + lx1m*ly1m*lz1m*lelt), stat=ierr) ! v2z
         if (ierr.ne.0) call exitti('alloc cb_gvolm$',ierr)
         cb_gvolm = 0

         allocate(cb_glog(
     $        1                        ! ifgeom
     $      + 1                        ! ifgmsh3
     $      + 1                        ! ifvcor
     $      + 1                        ! ifsurt
     $      + 1                        ! ifmelt
     $      + 1                        ! ifwcno
     $      + lelt                     ! ifrzer
     $      + 2*ldim*lelv              ! ifqinp
     $      + 2*ldim*lelv              ! ifeppm
     $      + 2                        ! iflmsf
     $      + 2                        ! iflmse
     $      + 2                        ! iflmsc
     $      + 2*ldim*lelt*2            ! ifmsfc
     $      + 12*lelt*2                ! ifmseg
     $      + 8*lelt*2                 ! ifmscr
     $      + 8*lelt                   ! ifnskp
     $      + 1), stat=ierr)           ! ifbcor
         if (ierr.ne.0) call exitti('alloc cb_glog$',ierr)
         cb_glog = .false.

         allocate(cb_cbbid(
     $        6*lelv                   ! boundaryID
     $      + 6*lelt), stat=ierr)      ! boundaryIDt
         if (ierr.ne.0) call exitti('alloc cb_cbbid$',ierr)
         cb_cbbid = 0

c        Group 1: /gxyz/
         ioff = 1
         ltot = lx1*ly1*lz1*lelt
         xm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gxyz(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         ym1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gxyz(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         zm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gxyz(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         ltot2 = lx2*ly2*lz2*lelv
         xm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_gxyz(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         ym2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_gxyz(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         zm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_gxyz(ioff : ioff + ltot2 - 1)

c        Group 2: /giso1/
         ioff = 1
         rxm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         sxm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         txm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         rym1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         sym1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         tym1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         rzm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         szm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         tzm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         jacm1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_giso1(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         jacmi(1:lx1*ly1*lz1,1:lelt) => cb_giso1(ioff : ioff + ltot - 1)

c        Group 3: /giso2/
         ioff = 1
         rxm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         sxm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         txm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         rym2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         sym2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         tym2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         rzm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         szm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         tzm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         jacm2(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_giso2(ioff : ioff + ltot2 - 1)

c        Group 4: /gisod/
         rx(1:lxd*lyd*lzd,1:ldim*ldim,1:lelv) =>
     $         cb_gisod(1 : lxd*lyd*lzd*ldim*ldim*lelv)

c        Group 5: /gmfact/
         ioff = 1
         g1m1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gmfact(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         g2m1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gmfact(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         g3m1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gmfact(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         g4m1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gmfact(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         g5m1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gmfact(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         g6m1(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_gmfact(ioff : ioff + ltot - 1)

c        Group 6: /gsurf/
         ioff = 1
         unr(1:lx1*lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         uns(1:lx1*lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         unt(1:lx1*lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         unx(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         uny(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         unz(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         t1x(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         t1y(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         t1z(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         t2x(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         t2y(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         t2z(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         area(1:lx1,1:lz1,1:6,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*6*lelt - 1)
         ioff = ioff + lx1*lz1*6*lelt
         etalph(1:lx1*lz1,1:2*ldim,1:lelt) =>
     $         cb_gsurf(ioff : ioff + lx1*lz1*2*ldim*lelt - 1)
         ioff = ioff + lx1*lz1*2*ldim*lelt
         dlam => cb_gsurf(ioff)

c        Group 7: /gvolm/
         ioff = 1
         vnx(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         vny(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         vnz(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         v1x(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         v1y(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         v1z(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         v2x(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         v2y(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)
         ioff = ioff + lx1m*ly1m*lz1m*lelt
         v2z(1:lx1m,1:ly1m,1:lz1m,1:lelt) =>
     $         cb_gvolm(ioff : ioff + lx1m*ly1m*lz1m*lelt - 1)

c        Group 8: /glog/
         ioff = 1
         ifgeom => cb_glog(ioff)
         ioff = ioff + 1
         ifgmsh3 => cb_glog(ioff)
         ioff = ioff + 1
         ifvcor => cb_glog(ioff)
         ioff = ioff + 1
         ifsurt => cb_glog(ioff)
         ioff = ioff + 1
         ifmelt => cb_glog(ioff)
         ioff = ioff + 1
         ifwcno => cb_glog(ioff)
         ioff = ioff + 1
         ifrzer(1:lelt) => cb_glog(ioff : ioff + lelt - 1)
         ioff = ioff + lelt
         ifqinp(1:2*ldim,1:lelv) =>
     $         cb_glog(ioff : ioff + 2*ldim*lelv - 1)
         ioff = ioff + 2*ldim*lelv
         ifeppm(1:2*ldim,1:lelv) =>
     $         cb_glog(ioff : ioff + 2*ldim*lelv - 1)
         ioff = ioff + 2*ldim*lelv
         iflmsf(0:1) => cb_glog(ioff : ioff + 2 - 1)
         ioff = ioff + 2
         iflmse(0:1) => cb_glog(ioff : ioff + 2 - 1)
         ioff = ioff + 2
         iflmsc(0:1) => cb_glog(ioff : ioff + 2 - 1)
         ioff = ioff + 2
         ifmsfc(1:2*ldim,1:lelt,0:1) =>
     $         cb_glog(ioff : ioff + 2*ldim*lelt*2 - 1)
         ioff = ioff + 2*ldim*lelt*2
         ifmseg(1:12,1:lelt,0:1) =>
     $         cb_glog(ioff : ioff + 12*lelt*2 - 1)
         ioff = ioff + 12*lelt*2
         ifmscr(1:8,1:lelt,0:1) =>
     $         cb_glog(ioff : ioff + 8*lelt*2 - 1)
         ioff = ioff + 8*lelt*2
         ifnskp(1:8,1:lelt) =>
     $         cb_glog(ioff : ioff + 8*lelt - 1)
         ioff = ioff + 8*lelt
         ifbcor => cb_glog(ioff)

c        Group 9: /cbbid/
         ioff = 1
         boundaryID(1:6,1:lelv) =>
     $         cb_cbbid(ioff : ioff + 6*lelv - 1)
         ioff = ioff + 6*lelv
         boundaryIDt(1:6,1:lelt) =>
     $         cb_cbbid(ioff : ioff + 6*lelt - 1)

      end subroutine init
      end module geom_mod
