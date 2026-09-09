      module gmres_mod

      implicit none
c
c     GMRES variables
c
c     w is a work vector
c     c and s store the Givens rotations
c     V stores the orthogonal Krylov subspace basis
c          -1
c     Z = M   V
c
      real, allocatable, target :: cb_gmres(:)
      real, allocatable, target :: cb_gmre1(:)
      real, allocatable, target :: cb_gmre2(:)
      real, allocatable, target :: cb_spltprec(:)

      real, pointer :: x_gmres(:), r_gmres(:)
      real, pointer :: w_gmres(:), h_gmres(:,:)
      real, pointer :: gamma_gmres(:), c_gmres(:), s_gmres(:)

      real, pointer :: v_gmres(:,:)

      real, pointer :: z_gmres(:,:)

      real, pointer :: ml_gmres(:), mu_gmres(:)

      contains

      subroutine init
         use size_mod
         implicit none
         integer ltot2

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_gmres(
     $        lx2*ly2*lz2*lelv       ! x_gmres
     $      + lx2*ly2*lz2*lelv       ! r_gmres
     $      + lx2*ly2*lz2*lelv       ! w_gmres
     $      + lgmres*lgmres          ! h_gmres
     $      + (lgmres+1)             ! gamma_gmres
     $      + lgmres                 ! c_gmres
     $      + lgmres), stat=ierr)    ! s_gmres
         if (ierr.ne.0) call exitti('alloc cb_gmres$',ierr)
         cb_gmres = 0

         allocate(cb_gmre1(lx2*ly2*lz2*lelv*lgmres), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_gmre1$',ierr)
         cb_gmre1 = 0
                                                       ! v_gmres

         allocate(cb_gmre2(lx2*ly2*lz2*lelv*lgmres), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_gmre2$',ierr)
         cb_gmre2 = 0
                                                       ! z_gmres

         allocate(cb_spltprec(
     $        lx2*ly2*lz2*lelv       ! ml_gmres
     $      + lx2*ly2*lz2*lelv), stat=ierr) ! mu_gmres
         if (ierr.ne.0) call exitti('alloc cb_spltprec$',ierr)
         cb_spltprec = 0

c        Group 1: /gmres/
         ioff = 1
         ltot2 = lx2*ly2*lz2*lelv
         x_gmres(1:lx2*ly2*lz2*lelv) =>
     $   cb_gmres(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         r_gmres(1:lx2*ly2*lz2*lelv) =>
     $   cb_gmres(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         w_gmres(1:lx2*ly2*lz2*lelv) =>
     $   cb_gmres(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         h_gmres(1:lgmres,1:lgmres) =>
     $         cb_gmres(ioff : ioff + lgmres*lgmres - 1)
         ioff = ioff + lgmres*lgmres
         gamma_gmres(1:lgmres+1) =>
     $         cb_gmres(ioff : ioff + (lgmres+1) - 1)
         ioff = ioff + (lgmres+1)
         c_gmres(1:lgmres) => cb_gmres(ioff : ioff + lgmres - 1)
         ioff = ioff + lgmres
         s_gmres(1:lgmres) => cb_gmres(ioff : ioff + lgmres - 1)

c        Group 2: /gmre1/
         v_gmres(1:lx2*ly2*lz2*lelv,1:lgmres) =>
     $         cb_gmre1(1 : lx2*ly2*lz2*lelv*lgmres)

c        Group 3: /gmre2/
         z_gmres(1:lx2*ly2*lz2*lelv,1:lgmres) =>
     $         cb_gmre2(1 : lx2*ly2*lz2*lelv*lgmres)

c        Group 4: /spltprec/
         ioff = 1
         ml_gmres(1:lx2*ly2*lz2*lelv) =>
     $   cb_spltprec(ioff : ioff + ltot2 - 1)
         ioff = ioff + ltot2
         mu_gmres(1:lx2*ly2*lz2*lelv) =>
     $   cb_spltprec(ioff : ioff + ltot2 - 1)

      end subroutine init
      end module gmres_mod
