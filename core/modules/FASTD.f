      module fastd_mod

      implicit none

      real, allocatable, target :: cb_fastd(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        fasts.f (local_solves_fdm) / navier6.f (set_overlap):
c        df(lx1*ly1*lz1,levb), sr,ss,st(lxx*2,levb); levb=lelv+lbelv,
c        lxx=lx1*lx1.  Computed once in navier6.f's gen_fast call,
c        read on every solve in fasts.f's local_solves_fdm.
         integer sz
         integer levb, lxx

         levb = lelv + lbelv
         lxx  = lx1*lx1

         sz = levb*(lx1*ly1*lz1 + 3*2*lxx)

         n = sz

         allocate(cb_fastd(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_fastd$',ierr)
         cb_fastd = 0

      end subroutine init
      end module fastd_mod
