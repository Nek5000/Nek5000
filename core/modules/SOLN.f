      module soln_mod

      implicit none

      real, allocatable, target    :: cb_bqcb(:)
      real, allocatable, target    :: cb_vptsol(:)
      real, allocatable, target    :: cb_vptsolm(:)
      real, allocatable, target    :: cb_expvis(:)
      real, allocatable, target    :: cb_cbm2(:)
      real, allocatable, target    :: cb_diverg(:)
      real, allocatable, target    :: cb_p0therm(:)
      real, allocatable, target    :: cb_vptmsk(:)
      real, allocatable, target    :: cb_pvptsl(:)
      integer, allocatable, target :: cb_ppointr(:)

      real, pointer :: bq(:,:,:,:,:), adq(:,:,:,:,:)
      real, pointer :: vxlag(:,:,:,:,:), vylag(:,:,:,:,:)
     $               , vzlag(:,:,:,:,:)
      real, pointer :: tlag(:,:,:,:,:,:), vgradt1(:,:,:,:,:)
     $               , vgradt2(:,:,:,:,:)
      real, pointer :: abx1(:,:,:,:), aby1(:,:,:,:), abz1(:,:,:,:)
      real, pointer :: abx2(:,:,:,:), aby2(:,:,:,:), abz2(:,:,:,:)
      real, pointer :: vdiff_e(:,:,:,:), vx(:,:,:,:)
     $               , vy(:,:,:,:), vz(:,:,:,:)
      real, pointer :: vx_e(:,:,:,:), vy_e(:,:,:,:), vz_e(:,:,:,:)
      real, pointer :: t(:,:,:,:,:), vtrans(:,:,:,:,:)
     $               , vdiff(:,:,:,:,:)
      real, pointer :: bfx(:,:,:,:), bfy(:,:,:,:), bfz(:,:,:,:)
     $               , cflf(:,:,:,:)
      real, pointer :: bmnv(:,:), bmass(:,:), bdivw(:,:), c_vx(:,:)
     $               , fw(:,:)
      real, pointer :: bx(:,:,:,:), by(:,:,:,:), bz(:,:,:,:)
     $               , pm(:,:,:,:)
      real, pointer :: bmx(:,:,:,:), bmy(:,:,:,:), bmz(:,:,:,:)
      real, pointer :: bbx1(:,:,:,:), bby1(:,:,:,:), bbz1(:,:,:,:)
      real, pointer :: bbx2(:,:,:,:), bby2(:,:,:,:), bbz2(:,:,:,:)
      real, pointer :: bxlag(:,:), bylag(:,:), bzlag(:,:), pmlag(:,:)
      real, pointer :: nu_star
      real, pointer :: pr(:,:,:,:), prlag(:,:,:,:,:)
      real, pointer :: qtl(:,:,:,:), usrdiv(:,:,:,:)
      real, pointer :: p0th, dp0thdt, gamma0, p0thn, p0thlag(:)
      real, pointer :: v1mask(:,:,:,:), v2mask(:,:,:,:), v3mask(:,:,:,:)
      real, pointer :: pmask(:,:,:,:), tmask(:,:,:,:,:), omask(:,:,:,:)
      real, pointer :: vmult(:,:,:,:), tmult(:,:,:,:,:)
      real, pointer :: b1mask(:,:,:,:), b2mask(:,:,:,:), b3mask(:,:,:,:)
     $               , bpmask(:,:,:,:)
      real, pointer :: vxp(:,:), vyp(:,:), vzp(:,:), prp(:,:)
     $               , tp(:,:,:)
      real, pointer :: bqp(:,:,:), adqp(:,:,:), bfxp(:,:), bfyp(:,:)
     $               , bfzp(:,:)
      real, pointer :: vxlagp(:,:,:), vylagp(:,:,:), vzlagp(:,:,:)
     $               , prlagp(:,:,:)
      real, pointer :: tlagp(:,:,:,:), exx1p(:,:), exy1p(:,:)
     $               , exz1p(:,:)
      real, pointer :: exx2p(:,:), exy2p(:,:), exz2p(:,:)
     $               , vgradt1p(:,:,:), vgradt2p(:,:,:)
      real, pointer :: adqp_dummy(:,:,:)

      integer, pointer :: jp


      contains

      subroutine init
         use size_mod
         implicit none
         integer ltot, ltotv, ltot2

         integer lorder2, ierr, ntotvsol, ntotvslp, ioff

         lorder2 = max(1, lorder-2)

c        --- allocate backing arrays ---

         allocate(cb_bqcb(lx1*ly1*lz1*lelt*ldimt * 2), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_bqcb$',ierr)
         cb_bqcb = 0

         ntotvsol =
     $      (lx1*ly1*lz1*lelv * 6)
     $    + (lx1*ly1*lz1*lelt*(lorder-1)*ldimt)
     $    + (lx1*ly1*lz1*lelt*ldimt * 2)
     $    + (lx1*ly1*lz1*lelv * 6)
     $    + (lx1*ly1*lz1*lelt)
     $    + (lx1*ly1*lz1*lelv * 3)
     $    + (lx1*ly1*lz1*lelt*ldimt)
     $    + (lx1*ly1*lz1*lelt*ldimt1 * 2)
     $    + (lx1*ly1*lz1*lelv * 4)
     $    + (lxd*lyd*lzd*lelv*ldim*(lorder+1))
     $    + (2*ldim*lelt)
     $    + (lx1*ly1*lz1*lelv*ldim*(lorder+1) * 3)
     $    + (lx1*ly1*lz1*lelv * 3)
         allocate(cb_vptsol(ntotvsol), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_vptsol$',ierr)
         cb_vptsol = 0

         allocate(cb_vptsolm(
     $      (lbx1*lby1*lbz1*lbelv * 3)
     $    + (lbx2*lby2*lbz2*lbelv)
     $    + (lbx1*lby1*lbz1*lbelv * 9)
     $    + (lbx1*lby1*lbz1*lbelv*(lorder-1) * 3)
     $    + (lbx2*lby2*lbz2*lbelv*lorder2)), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_vptsolm$',ierr)
         cb_vptsolm = 0

         allocate(cb_expvis(1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_expvis$',ierr)
         cb_expvis = 0
         allocate(cb_cbm2(
     $      (lx2*ly2*lz2*lelv)
     $    + (lx2*ly2*lz2*lelv*lorder2)), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_cbm2$',ierr)
         cb_cbm2 = 0
         allocate(cb_diverg(lx2*ly2*lz2*lelt * 2), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_diverg$',ierr)
         cb_diverg = 0
         allocate(cb_p0therm(6), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_p0therm$',ierr)
         cb_p0therm = 0

         allocate(cb_vptmsk(
     $      (lx1*ly1*lz1*lelv * 4)
     $    + (lx1*ly1*lz1*lelt*ldimt)
     $    + (lx1*ly1*lz1*lelt)
     $    + (lx1*ly1*lz1*lelv)
     $    + (lx1*ly1*lz1*lelt*ldimt)
     $    + (lbx1*lby1*lbz1*lbelv * 4)), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_vptmsk$',ierr)
         cb_vptmsk = 0

         ntotvslp =
     $      (lpx1*lpy1*lpz1*lpelv*lpert * 3)
     $    + (lpx2*lpy2*lpz2*lpelv*lpert)
     $    + (lpx1*lpy1*lpz1*lpelt*ldimt*lpert * 3)
     $    + (lpx1*lpy1*lpz1*lpelv*lpert * 3)
     $    + (lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert * 3)
     $    + (lpx2*lpy2*lpz2*lpelv*lorder2*lpert)
     $    + (lpx1*lpy1*lpz1*lpelt*ldimt*(lorder-1)*lpert)
     $    + (lpx1*lpy1*lpz1*lpelv*lpert * 6)
     $    + (lpx1*lpy1*lpz1*lpelt*ldimt*lpert * 2)
         allocate(cb_pvptsl(ntotvslp), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_pvptsl$',ierr)
         cb_pvptsl = 0

         allocate(cb_ppointr(1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_ppointr$',ierr)
         cb_ppointr = 0

c        Group 1: /bqcb/
         ioff = 1
         bq(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_bqcb(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt
         adq(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_bqcb(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)

c        Group 2: /vptsol/
         ioff = 1
         vxlag(1:lx1,1:ly1,1:lz1,1:lelv,1:2) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelv*2 - 1)
         ioff = ioff + lx1*ly1*lz1*lelv*2
         vylag(1:lx1,1:ly1,1:lz1,1:lelv,1:2) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelv*2 - 1)
         ioff = ioff + lx1*ly1*lz1*lelv*2
         vzlag(1:lx1,1:ly1,1:lz1,1:lelv,1:2) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelv*2 - 1)
         ioff = ioff + lx1*ly1*lz1*lelv*2
         tlag(1:lx1,1:ly1,1:lz1,1:lelt,1:lorder-1,1:ldimt) =>
     $         cb_vptsol(ioff : ioff
     $         + lx1*ly1*lz1*lelt*(lorder-1)*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*(lorder-1)*ldimt
         vgradt1(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt
         vgradt2(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt
         ltotv = lx1*ly1*lz1*lelv
         abx1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         aby1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         abz1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         abx2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         aby2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         abz2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         ltot = lx1*ly1*lz1*lelt
         vdiff_e(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_vptsol(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         vx(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vy(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vz(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         t(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt
         vtrans(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt1) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelt*ldimt1 - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt1
         vdiff(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt1) =>
     $         cb_vptsol(ioff : ioff + lx1*ly1*lz1*lelt*ldimt1 - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt1
         bfx(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         bfy(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         bfz(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         cflf(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         c_vx(1:lxd*lyd*lzd*lelv*ldim,1:lorder+1) =>
     $         cb_vptsol(ioff : ioff
     $         + lxd*lyd*lzd*lelv*ldim*(lorder+1) - 1)
         ioff = ioff + lxd*lyd*lzd*lelv*ldim*(lorder+1)
         fw(1:2*ldim,1:lelt) =>
     $         cb_vptsol(ioff : ioff + 2*ldim*lelt - 1)
         ioff = ioff + 2*ldim*lelt
         bmnv(1:lx1*ly1*lz1*lelv*ldim,1:lorder+1) =>
     $         cb_vptsol(ioff : ioff
     $         + lx1*ly1*lz1*lelv*ldim*(lorder+1) - 1)
         ioff = ioff + lx1*ly1*lz1*lelv*ldim*(lorder+1)
         bmass(1:lx1*ly1*lz1*lelv*ldim,1:lorder+1) =>
     $         cb_vptsol(ioff : ioff
     $         + lx1*ly1*lz1*lelv*ldim*(lorder+1) - 1)
         ioff = ioff + lx1*ly1*lz1*lelv*ldim*(lorder+1)
         bdivw(1:lx1*ly1*lz1*lelv*ldim,1:lorder+1) =>
     $         cb_vptsol(ioff : ioff
     $         + lx1*ly1*lz1*lelv*ldim*(lorder+1) - 1)
         ioff = ioff + lx1*ly1*lz1*lelv*ldim*(lorder+1)
         vx_e(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vy_e(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vz_e(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptsol(ioff : ioff + ltotv - 1)

c        Group 3: /vptsolm/
         ioff = 1
         bx(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         by(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bz(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         pm(1:lbx2,1:lby2,1:lbz2,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx2*lby2*lbz2*lbelv - 1)
         ioff = ioff + lbx2*lby2*lbz2*lbelv
         bmx(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bmy(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bmz(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bbx1(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bby1(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bbz1(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bbx2(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bby2(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bbz2(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptsolm(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bxlag(1:lbx1*lby1*lbz1*lbelv,1:lorder-1) =>
     $         cb_vptsolm(ioff : ioff
     $         + lbx1*lby1*lbz1*lbelv*(lorder-1) - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv*(lorder-1)
         bylag(1:lbx1*lby1*lbz1*lbelv,1:lorder-1) =>
     $         cb_vptsolm(ioff : ioff
     $         + lbx1*lby1*lbz1*lbelv*(lorder-1) - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv*(lorder-1)
         bzlag(1:lbx1*lby1*lbz1*lbelv,1:lorder-1) =>
     $         cb_vptsolm(ioff : ioff
     $         + lbx1*lby1*lbz1*lbelv*(lorder-1) - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv*(lorder-1)
         pmlag(1:lbx2*lby2*lbz2*lbelv,1:lorder2) =>
     $         cb_vptsolm(ioff : ioff
     $         + lbx2*lby2*lbz2*lbelv*lorder2 - 1)

c        Group 4: /expvis/, /cbm2/, /diverg/, /p0therm/
         nu_star => cb_expvis(1)
         ioff = 1
         ltot2 = lx2*ly2*lz2*lelv
         pr(1:lx2,1:ly2,1:lz2,1:lelv) =>
     $   cb_cbm2(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         prlag(1:lx2,1:ly2,1:lz2,1:lelv,1:lorder2) =>
     $         cb_cbm2(ioff : ioff + lx2*ly2*lz2*lelv*lorder2 - 1)
         ioff = 1
         qtl(1:lx2,1:ly2,1:lz2,1:lelt) =>
     $         cb_diverg(ioff : ioff + lx2*ly2*lz2*lelt - 1)
         ioff = ioff + lx2*ly2*lz2*lelt
         usrdiv(1:lx2,1:ly2,1:lz2,1:lelt) =>
     $         cb_diverg(ioff : ioff + lx2*ly2*lz2*lelt - 1)
         p0th    => cb_p0therm(1)
         dp0thdt => cb_p0therm(2)
         gamma0  => cb_p0therm(3)
         p0thn   => cb_p0therm(4)
         p0thlag => cb_p0therm(5 : 6)

c        Group 8: /vptmsk/
         ioff = 1
         v1mask(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptmsk(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         v2mask(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptmsk(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         v3mask(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptmsk(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         pmask(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptmsk(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         tmask(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_vptmsk(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt
         omask(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_vptmsk(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         vmult(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_vptmsk(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         tmult(1:lx1,1:ly1,1:lz1,1:lelt,1:ldimt) =>
     $         cb_vptmsk(ioff : ioff + lx1*ly1*lz1*lelt*ldimt - 1)
         ioff = ioff + lx1*ly1*lz1*lelt*ldimt
         b1mask(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptmsk(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         b2mask(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptmsk(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         b3mask(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptmsk(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)
         ioff = ioff + lbx1*lby1*lbz1*lbelv
         bpmask(1:lbx1,1:lby1,1:lbz1,1:lbelv) =>
     $         cb_vptmsk(ioff : ioff + lbx1*lby1*lbz1*lbelv - 1)

c        Group 9: /pvptsl/
         ioff = 1
         vxp(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         vyp(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         vzp(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         prp(1:lpx2*lpy2*lpz2*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx2*lpy2*lpz2*lpelv*lpert - 1)
         ioff = ioff + lpx2*lpy2*lpz2*lpelv*lpert
         tp(1:lpx1*lpy1*lpz1*lpelt,1:ldimt,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelt*ldimt*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelt*ldimt*lpert
         bqp(1:lpx1*lpy1*lpz1*lpelt,1:ldimt,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelt*ldimt*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelt*ldimt*lpert
         bfxp(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         bfyp(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         bfzp(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         vxlagp(1:lpx1*lpy1*lpz1*lpelv,1:lorder-1,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert
         vylagp(1:lpx1*lpy1*lpz1*lpelv,1:lorder-1,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert
         vzlagp(1:lpx1*lpy1*lpz1*lpelv,1:lorder-1,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*(lorder-1)*lpert
         prlagp(1:lpx2*lpy2*lpz2*lpelv,1:lorder2,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx2*lpy2*lpz2*lpelv*lorder2*lpert - 1)
         ioff = ioff + lpx2*lpy2*lpz2*lpelv*lorder2*lpert
         tlagp(1:lpx1*lpy1*lpz1*lpelt,1:ldimt,1:lorder-1,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelt*ldimt*(lorder-1)*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelt*ldimt*(lorder-1)*lpert
         exx1p(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         exy1p(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         exz1p(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         exx2p(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         exy2p(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         exz2p(1:lpx1*lpy1*lpz1*lpelv,1:lpert) =>
     $         cb_pvptsl(ioff : ioff + lpx1*lpy1*lpz1*lpelv*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelv*lpert
         vgradt1p(1:lpx1*lpy1*lpz1*lpelt,1:ldimt,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelt*ldimt*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelt*ldimt*lpert
         vgradt2p(1:lpx1*lpy1*lpz1*lpelt,1:ldimt,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelt*ldimt*lpert - 1)
         ioff = ioff + lpx1*lpy1*lpz1*lpelt*ldimt*lpert
         adqp(1:lpx1*lpy1*lpz1*lpelt,1:ldimt,1:lpert) =>
     $         cb_pvptsl(ioff : ioff
     $         + lpx1*lpy1*lpz1*lpelt*ldimt*lpert - 1)

c        Group 10: /ppointr/
         jp => cb_ppointr(1)

      end subroutine init
      end module soln_mod
