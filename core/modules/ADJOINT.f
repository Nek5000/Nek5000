      module adjoint_mod

      implicit none
c
c     Adjoint variables
c
      real,    allocatable, target :: cb_adj_real(:)
      integer, allocatable, target :: cb_adj_integer(:)
      logical, allocatable, target :: cb_adj_logical(:)
      real,    allocatable, target :: cb_dtgrad(:)
      real,    allocatable, target :: cb_gravity_adjoint(:)

      real, pointer :: vxadj(:,:,:,:), vyadj(:,:,:,:), vzadj(:,:,:,:)
      real, pointer :: vxadjold(:,:,:,:), vyadjold(:,:,:,:)
      real, pointer :: vzadjold(:,:,:,:)
      real, pointer :: endtime, adjtol
      real, pointer :: tpadj(:,:,:,:), tpadjold(:,:,:,:)
      real, pointer :: alpha_max

      integer, pointer :: npassadj, maxpassadj

      logical, pointer :: ifadj

      real, pointer :: dTdx(:), dTdy(:), dTdz(:)

      real, pointer :: g_adj(:)
      real, pointer :: beta_b(:)

      contains

      subroutine init
         use size_mod
         implicit none
         integer ltot, ltotv

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_adj_real(
     $        lx1*ly1*lz1*lelv         ! vxadj
     $      + lx1*ly1*lz1*lelv         ! vyadj
     $      + lx1*ly1*lz1*lelv         ! vzadj
     $      + lx1*ly1*lz1*lelv         ! vxadjold
     $      + lx1*ly1*lz1*lelv         ! vyadjold
     $      + lx1*ly1*lz1*lelv         ! vzadjold
     $      + 1                        ! endtime
     $      + 1                        ! adjtol
     $      + lx1*ly1*lz1*lelt         ! tpadj
     $      + lx1*ly1*lz1*lelt         ! tpadjold
     $      + 1), stat=ierr)           ! alpha_max
         if (ierr.ne.0) call exitti('alloc cb_adj_real$',ierr)
         cb_adj_real = 0

         allocate(cb_adj_integer(2), stat=ierr) ! npassadj, maxpassadj
         if (ierr.ne.0) call exitti('alloc cb_adj_integer$',ierr)
         cb_adj_integer = 0

         allocate(cb_adj_logical(1), stat=ierr) ! ifadj
         if (ierr.ne.0) call exitti('alloc cb_adj_logical$',ierr)
         cb_adj_logical = .false.

         allocate(cb_dtgrad(
     $        lx1*ly1*lz1*lelt         ! dTdx
     $      + lx1*ly1*lz1*lelt         ! dTdy
     $      + lx1*ly1*lz1*lelt), stat=ierr) ! dTdz
         if (ierr.ne.0) call exitti('alloc cb_dtgrad$',ierr)
         cb_dtgrad = 0

         allocate(cb_gravity_adjoint(
     $        3                        ! g_adj
     $      + lx1*ly1*lz1*lelt), stat=ierr) ! beta_b
         if (ierr.ne.0) call exitti('alloc cb_gravity_adjoint$',ierr)
         cb_gravity_adjoint = 0

c        Group 1: /adj_real/
         ioff = 1
         ltotv = lx1*ly1*lz1*lelv
         vxadj(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_adj_real(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vyadj(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_adj_real(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vzadj(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_adj_real(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vxadjold(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_adj_real(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vyadjold(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_adj_real(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         vzadjold(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_adj_real(ioff : ioff + ltotv - 1)
         ioff = ioff + ltotv
         endtime => cb_adj_real(ioff)
         ioff = ioff + 1
         adjtol => cb_adj_real(ioff)
         ioff = ioff + 1
         ltot = lx1*ly1*lz1*lelt
         tpadj(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_adj_real(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         tpadjold(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $   cb_adj_real(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         alpha_max => cb_adj_real(ioff)

c        Group 2: /adj_integer/
         ioff = 1
         npassadj => cb_adj_integer(ioff)
         ioff = ioff + 1
         maxpassadj => cb_adj_integer(ioff)

c        Group 3: /adj_logical/
         ifadj => cb_adj_logical(1)

c        Group 4: /dTgrad/
         ioff = 1
         dTdx(1:lx1*ly1*lz1*lelt) => cb_dtgrad(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         dTdy(1:lx1*ly1*lz1*lelt) => cb_dtgrad(ioff : ioff + ltot - 1)
         ioff = ioff + ltot
         dTdz(1:lx1*ly1*lz1*lelt) => cb_dtgrad(ioff : ioff + ltot - 1)

c        Group 5: /gravity_adjoint/
         ioff = 1
         g_adj(1:3) => cb_gravity_adjoint(ioff : ioff + 3 - 1)
         ioff = ioff + 3
         beta_b(1:lx1*ly1*lz1*lelt) =>
     $   cb_gravity_adjoint(ioff : ioff + ltot - 1)

      end subroutine init
      end module adjoint_mod
