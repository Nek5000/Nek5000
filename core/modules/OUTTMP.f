      module outtmp_mod

      implicit none

      real, allocatable, target :: cb_outtmp(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        postpro.f (hpts: wrk(lx1*ly1*lz1*lelt,nfldm), nfldm=ldim+ldimt+1)
         integer sz_hpts
c        prepost.f (outpost2: w1,w2,w3(ltot1),wp(ltot2),wt(ltot1,ldimt))
         integer sz_outpost2
         integer ltot1, ltot2

         ltot1 = lx1*ly1*lz1*lelt
         ltot2 = lx2*ly2*lz2*lelv

         sz_hpts     = ltot1*(ldim+ldimt+1)
         sz_outpost2 = ltot1*(3+ldimt) + ltot2

         n = max(sz_hpts, sz_outpost2)

         allocate(cb_outtmp(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_outtmp$',ierr)
         cb_outtmp = 0

      end subroutine init
      end module outtmp_mod
