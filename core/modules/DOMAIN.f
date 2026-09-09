      module domain_mod
      include 'SIZE'
c
c     Arrays for overlapping Schwartz algorithm
c
      integer ltotd
      parameter (ltotd = lx1*ly1*lz1*lelt                     )
c
      integer ndom, n_o, nel_proc, gs_hnd_overlap
      common /ddptri/ ndom,n_o,nel_proc,gs_hnd_overlap

      integer, allocatable, target :: na(:), ma(:), nza(:)
c
c     These are the H1 coarse-grid arrays:
c
      integer lxc, lcr
      parameter(lxc=2)
      parameter(lcr=lxc**ldim)

      integer*8, allocatable, target :: se_to_gcrs(:,:)
      integer n_crs,m_crs,nx_crs,nxyz_c
      common /h1_crsi/ n_crs,m_crs, nx_crs, nxyz_c
c
      real             h1_basis(lx1*lxc), h1_basist(lxc*lx1)
      common /h1_crs/  h1_basis         , h1_basist

      real             l2_basis(lx2*lxc), l2_basist(lxc*lx2)
      equivalence     (h1_basis  , l2_basis  )
      equivalence     (h1_basist , l2_basist )

      contains

      subroutine init
         implicit none
         integer ierr

         allocate(na(lelt+1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc na$',ierr)
         na = 0
         allocate(ma(lelt+1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc ma$',ierr)
         ma = 0
         allocate(nza(lelt+1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc nza$',ierr)
         nza = 0
         allocate(se_to_gcrs(lcr,lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc se_to_gcrs$',ierr)
         se_to_gcrs = 0_8

      end subroutine init
      end module domain_mod
