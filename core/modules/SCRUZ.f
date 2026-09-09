      module scruz_mod

      implicit none

      real, allocatable, target :: cb_scruz(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscruz
c        rans_komg.f (comp_StOm: sij(*,6,lelv),oij(*,lelv,3));
c        also covers plan4.f (crespsp: sij(lx1*ly1*lz1,6,lelv))
         integer sz_sij
         integer sz_ur    ! prepost.f (mfo_outfld: ur1,ur2,ur3(lxo**3,lelt))
         integer sz_mid   ! postpro.f/navier5.f (mid(2*lt)/tmsk,tmlt,w1,w2)
         integer sz_b     ! mxm_std.f (mxmtest: b(4*lx1*ly1*lz1*lelt))
         integer sz_fix   ! FIXED-size views:
c                           comm_mpi.f gop_test times(2,500)=1000 real;
c                           postpro.f RE2 tools igr(1000), icurve(12,500)=6000 int.

         sz_sij = 9*lx1*ly1*lz1*lelv
         sz_ur  = 3*lxo*lxo*lxo*lelt
         sz_mid = 2*lx1*ly1*lz1*lelt
         sz_b   = 4*lx1*ly1*lz1*lelt
         sz_fix = 12*500     ! icurve(12,500) integer view (6000 int >= all others)

         nscruz = max(sz_sij, sz_ur, sz_mid, sz_b, sz_fix)

         allocate(cb_scruz(nscruz), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scruz$',ierr)
         cb_scruz = 0

      end subroutine init
      end module scruz_mod
