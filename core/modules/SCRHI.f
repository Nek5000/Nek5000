      module scrhi_mod

      implicit none

      real, allocatable, target :: cb_scrhi(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        h2inv(lx1,ly1,lz1,lelv) at 8 sites (multimesh.f, navier1.f,
c        perturb.f, eigsolv.f, induct.f, drive2.f, hsmg.f x2);
c        h2inv(lx1,ly1,lz1,lelt) at hsmg.f's h1mg_setup (dominates,
c        lelt>=lelv)
         integer sz

         sz = lx1*ly1*lz1*lelt

         n = sz

         allocate(cb_scrhi(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrhi$',ierr)
         cb_scrhi = 0

      end subroutine init
      end module scrhi_mod
