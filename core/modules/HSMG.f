      module hsmg_mod

      implicit none
c
c     Variables for the Hybrid Schwarz Multigrid
c
c     Allocate MHD memory only if lbx1==lx1
c
c     lmg_mhd    : 1 if MHD is true, 0 otherwise
c     lmgs       : max number of multigrid solvers
c     lmgn       : max number of multigrid levels
c     lmgx       : max number of mg index levels
c     lxm,lym,lzm: mgrid sizes
c     lmg_rwt    : restriction weight max size
c     lmg_fasts  : FDM S max size
c     lmg_fastd  : FDM D max size
c     lmg_swt    : schwarz weight max size
c     lmg_g      : metrics max size
c     lmg_solve  : solver r,e max size
c
c     set in init() below (depend on SIZE, so cannot be PARAMETERs here)
c
      integer lmg_mhd,lmgs,lmgn,lmgx,lxm,lym,lzm,lmg_rwt,lmg_fasts,
     $        lmg_fastd,lmg_swt,lmg_g,lmg_solve

      real, allocatable, target    :: cb_mghr(:)
      integer, allocatable, target :: cb_mghs(:)
      integer, allocatable, target :: cb_mgh1i(:)

      integer, pointer :: mg_lmax                 !number of multigrid levels
      integer, pointer :: mg_nx(:)                !level poly. order (for GLL pts)
      integer, pointer :: mg_ny(:), mg_nz(:)
      integer, pointer :: mg_nh(:), mg_nhz(:)      !number of 1d nodes
      integer, pointer :: mg_gsh_schwarz_handle(:,:) !dssum schwarz handles
      integer, pointer :: mg_gsh_handle(:,:)          !dssum handle
      integer, pointer :: mg_rstr_wt_index(:,:)
      integer, pointer :: mg_mask_index(:,:)
      integer, pointer :: mg_solve_index(:,:)
      integer, pointer :: mg_fast_s_index(:,:)
      integer, pointer :: mg_fast_d_index(:,:)
      integer, pointer :: mg_schwarz_wt_index(:,:)
      integer, pointer :: mg_g_index(:,:)
      integer, pointer :: mg_fld                  !active mg field

      real, pointer :: mg_jh(:,:)      !c-to-f interpolation matrices
      real, pointer :: mg_jht(:,:)     !transpose of mg_jh
      real, pointer :: mg_jhfc(:,:)    !f-to-c interpolation
      real, pointer :: mg_jhfct(:,:)   !transpose of mg_jhfc
      real, pointer :: mg_ah(:,:)      !A hat matrices
      real, pointer :: mg_bh(:,:)      !B hat matrices
      real, pointer :: mg_dh(:,:)      !D hat matrices
      real, pointer :: mg_dht(:,:)     !D hat transpose matrices
      real, pointer :: mg_zh(:,:)      !Nodal coordinates

      real, pointer :: mg_rstr_wt(:)   !restriction wt
      real, pointer :: mg_mask(:)      !b.c. mask
      real, pointer :: mg_fast_s(:)
      real, pointer :: mg_fast_d(:)
      real, pointer :: mg_schwarz_wt(:)
      real, pointer :: mg_solve_e(:)
      real, pointer :: mg_solve_r(:)

      real, pointer :: mg_h1(:)
      real, pointer :: mg_h2(:)
      real, pointer :: mg_b(:)
      real, pointer :: mg_g(:)         !metrics matrices

      real, pointer :: mg_work(:)      ! must be able to hold
      real, pointer :: mg_work2(:)     ! two lower level extended
      real, pointer :: mg_worke(:,:)   ! schwarz arrays

      integer, pointer :: mg_imask(:)  ! For h1mg, mask is a ptr

c     Specific to h1 multigrid:

      integer, pointer :: mg_h1_lmax
      integer, pointer :: mg_h1_n(:,:)
      integer, pointer :: p_mg_h1(:,:), p_mg_h2(:,:)
      integer, pointer :: p_mg_b(:,:), p_mg_g(:,:)
      integer, pointer :: p_mg_msk(:,:)

      contains

      subroutine init
         use size_mod
         use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
         implicit none

         integer ierr, ioff
         integer, pointer :: p_imask(:)

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

c        --- set mg sizing parameters (depend on SIZE) ---

         lmg_mhd   = 1-(lx1-lbx1)/(lx1-1) !1 if MHD is true, 0 otherwise
         lmgs      = 1 + lmg_mhd         ! max number of multigrid solvers
         lmgn      = 3                   ! max number of multigrid levels
         lmgx      = lmgn+1              ! max number of mg index levels
         lxm       = lx2+2               ! mgrid sizes
         lym       = lxm
         lzm       = lz2+2*(ldim-2)
         lmg_rwt   = 2*lxm*lzm           ! restriction weight max size
         lmg_fasts = 2*lxm*lxm           ! FDM S max size
         lmg_fastd = 2*lxm*lym*lzm       ! FDM D max size
         lmg_swt   = 2*lxm*lzm           ! schwarz weight max size
         lmg_g     = 2*lx2*ly2*lz2       ! metrics max size
         lmg_solve = 2*lxm*lym*lzm       ! solver r,e max size

c        --- allocate backing arrays ---

         allocate(cb_mghs(
     $        1                     ! mg_lmax
     $      + lmgn                  ! mg_nx
     $      + lmgn                  ! mg_ny
     $      + lmgn                  ! mg_nz
     $      + lmgn                  ! mg_nh
     $      + lmgn                  ! mg_nhz
     $      + lmgn*lmgs             ! mg_gsh_schwarz_handle
     $      + lmgn*lmgs             ! mg_gsh_handle
     $      + lmgx*(lmgs+1)         ! mg_rstr_wt_index
     $      + lmgx*(lmgs+1)         ! mg_mask_index
     $      + lmgx*(lmgs+1)         ! mg_solve_index
     $      + lmgx*(lmgs+1)         ! mg_fast_s_index
     $      + lmgx*(lmgs+1)         ! mg_fast_d_index
     $      + lmgx*(lmgs+1)         ! mg_schwarz_wt_index
     $      + lmgx*(lmgs+1)         ! mg_g_index
     $      + 1), stat=ierr)        ! mg_fld
         if (ierr.ne.0) call exitti('alloc cb_mghs$',ierr)
         cb_mghs = 0

         allocate(cb_mghr(
     $        lxm*lxm*lmgn                   ! mg_jh
     $      + lxm*lxm*lmgn                   ! mg_jht
     $      + lxm*lxm*lmgn                   ! mg_ah
     $      + lxm*lmgn                       ! mg_bh
     $      + lxm*lxm*lmgn                   ! mg_dh
     $      + lxm*lxm*lmgn                   ! mg_dht
     $      + lxm*lmgn                       ! mg_zh
     $      + lxm*lxm*lmgn                   ! mg_jhfc
     $      + lxm*lxm*lmgn                   ! mg_jhfct
     $      + lmgs*lmg_rwt*2*ldim*lelt       ! mg_rstr_wt
     $      + lmgs*lmg_rwt*4*ldim*lelt       ! mg_mask
     $      + lmgs*lmg_fasts*2*ldim*lelt     ! mg_fast_s
     $      + lmgs*lmg_fastd*lelt            ! mg_fast_d
     $      + lmgs*lmg_swt*4*ldim*lelt       ! mg_schwarz_wt
     $      + lmg_solve*lelt                 ! mg_solve_e
     $      + lmg_solve*lelt                 ! mg_solve_r
     $      + lmg_g*lelt                     ! mg_h1
     $      + lmg_g*lelt                     ! mg_h2
     $      + lmg_g*lelt                     ! mg_b
     $      + lmg_g*((ldim-1)*3)*lelt        ! mg_g
     $      + 2*lxm*lym*lzm*lelt             ! mg_work
     $      + lxm*lym*lzm*lelt               ! mg_work2
     $      + lxm*lym*lzm*6), stat=ierr)     ! mg_worke
         if (ierr.ne.0) call exitti('alloc cb_mghr$',ierr)
         cb_mghr = 0

         allocate(cb_mgh1i(
     $        1                     ! mg_h1_lmax
     $      + lmgx*ldimt1           ! mg_h1_n
     $      + lmgx*ldimt1           ! p_mg_h1
     $      + lmgx*ldimt1           ! p_mg_h2
     $      + lmgx*ldimt1           ! p_mg_b
     $      + lmgx*ldimt1           ! p_mg_g
     $      + lmgx*ldimt1), stat=ierr) ! p_mg_msk
         if (ierr.ne.0) call exitti('alloc cb_mgh1i$',ierr)
         cb_mgh1i = 0

c        Group 1: /mghs/
         ioff = 1
         mg_lmax => cb_mghs(ioff)
         ioff = ioff + 1
         mg_nx(1:lmgn) => cb_mghs(ioff : ioff + lmgn - 1)
         ioff = ioff + lmgn
         mg_ny(1:lmgn) => cb_mghs(ioff : ioff + lmgn - 1)
         ioff = ioff + lmgn
         mg_nz(1:lmgn) => cb_mghs(ioff : ioff + lmgn - 1)
         ioff = ioff + lmgn
         mg_nh(1:lmgn) => cb_mghs(ioff : ioff + lmgn - 1)
         ioff = ioff + lmgn
         mg_nhz(1:lmgn) => cb_mghs(ioff : ioff + lmgn - 1)
         ioff = ioff + lmgn
         mg_gsh_schwarz_handle(1:lmgn,1:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgn*lmgs - 1)
         ioff = ioff + lmgn*lmgs
         mg_gsh_handle(1:lmgn,1:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgn*lmgs - 1)
         ioff = ioff + lmgn*lmgs
         mg_rstr_wt_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_mask_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_solve_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_fast_s_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_fast_d_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_schwarz_wt_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_g_index(1:lmgx,0:lmgs) =>
     $         cb_mghs(ioff : ioff + lmgx*(lmgs+1) - 1)
         ioff = ioff + lmgx*(lmgs+1)
         mg_fld => cb_mghs(ioff)

c        Group 2: /mghr/
         ioff = 1
         mg_jh(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_jht(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_ah(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_bh(1:lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lmgn - 1)
         ioff = ioff + lxm*lmgn
         mg_dh(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_dht(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_zh(1:lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lmgn - 1)
         ioff = ioff + lxm*lmgn
         mg_jhfc(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_jhfct(1:lxm*lxm,1:lmgn) =>
     $         cb_mghr(ioff : ioff + lxm*lxm*lmgn - 1)
         ioff = ioff + lxm*lxm*lmgn
         mg_rstr_wt(0 : lmgs*lmg_rwt*2*ldim*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmgs*lmg_rwt*2*ldim*lelt - 1)
         ioff = ioff + lmgs*lmg_rwt*2*ldim*lelt
         mg_mask(0 : lmgs*lmg_rwt*4*ldim*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmgs*lmg_rwt*4*ldim*lelt - 1)
         ioff = ioff + lmgs*lmg_rwt*4*ldim*lelt
         mg_fast_s(0 : lmgs*lmg_fasts*2*ldim*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmgs*lmg_fasts*2*ldim*lelt - 1)
         ioff = ioff + lmgs*lmg_fasts*2*ldim*lelt
         mg_fast_d(0 : lmgs*lmg_fastd*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmgs*lmg_fastd*lelt - 1)
         ioff = ioff + lmgs*lmg_fastd*lelt
         mg_schwarz_wt(0 : lmgs*lmg_swt*4*ldim*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmgs*lmg_swt*4*ldim*lelt - 1)
         ioff = ioff + lmgs*lmg_swt*4*ldim*lelt
         mg_solve_e(0 : lmg_solve*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmg_solve*lelt - 1)
         ioff = ioff + lmg_solve*lelt
         mg_solve_r(0 : lmg_solve*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmg_solve*lelt - 1)
         ioff = ioff + lmg_solve*lelt
         mg_h1(0 : lmg_g*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmg_g*lelt - 1)
         ioff = ioff + lmg_g*lelt
         mg_h2(0 : lmg_g*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmg_g*lelt - 1)
         ioff = ioff + lmg_g*lelt
         mg_b(0 : lmg_g*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmg_g*lelt - 1)
         ioff = ioff + lmg_g*lelt
         mg_g(0 : lmg_g*((ldim-1)*3)*lelt - 1) =>
     $         cb_mghr(ioff : ioff + lmg_g*((ldim-1)*3)*lelt - 1)
         ioff = ioff + lmg_g*((ldim-1)*3)*lelt
         mg_work(1 : 2*lxm*lym*lzm*lelt) =>
     $         cb_mghr(ioff : ioff + 2*lxm*lym*lzm*lelt - 1)
         ioff = ioff + 2*lxm*lym*lzm*lelt
         mg_work2(1 : lxm*lym*lzm*lelt) =>
     $         cb_mghr(ioff : ioff + lxm*lym*lzm*lelt - 1)
         ioff = ioff + lxm*lym*lzm*lelt
         mg_worke(1:lxm*lym*lzm,1:6) =>
     $         cb_mghr(ioff : ioff + lxm*lym*lzm*6 - 1)

c        mg_imask: integer alias of mg_mask storage (former EQUIVALENCE)
         call c_f_pointer(c_loc(mg_mask(0)), p_imask,
     $         [lmgs*lmg_rwt*4*ldim*lelt])
         mg_imask(0 : lmgs*lmg_rwt*4*ldim*lelt - 1) => p_imask

c        Group 3: /mgh1i/
         ioff = 1
         mg_h1_lmax => cb_mgh1i(ioff)
         ioff = ioff + 1
         mg_h1_n(1:lmgx,1:ldimt1) =>
     $         cb_mgh1i(ioff : ioff + lmgx*ldimt1 - 1)
         ioff = ioff + lmgx*ldimt1
         p_mg_h1(1:lmgx,1:ldimt1) =>
     $         cb_mgh1i(ioff : ioff + lmgx*ldimt1 - 1)
         ioff = ioff + lmgx*ldimt1
         p_mg_h2(1:lmgx,1:ldimt1) =>
     $         cb_mgh1i(ioff : ioff + lmgx*ldimt1 - 1)
         ioff = ioff + lmgx*ldimt1
         p_mg_b(1:lmgx,1:ldimt1) =>
     $         cb_mgh1i(ioff : ioff + lmgx*ldimt1 - 1)
         ioff = ioff + lmgx*ldimt1
         p_mg_g(1:lmgx,1:ldimt1) =>
     $         cb_mgh1i(ioff : ioff + lmgx*ldimt1 - 1)
         ioff = ioff + lmgx*ldimt1
         p_mg_msk(1:lmgx,1:ldimt1) =>
     $         cb_mgh1i(ioff : ioff + lmgx*ldimt1 - 1)

      end subroutine init
      end module hsmg_mod
