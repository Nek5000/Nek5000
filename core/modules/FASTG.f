      module fastg_mod

      implicit none

      real, allocatable, target :: cb_fastg(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        gmres.f (gen_fast_g/fastdm1_g): sr,ss,st(lxss,2,lelv),
c        df(lxs*lys*lzs,lelv); lxss=lxs*lxs, lxs=lys=lzs=1
         integer sz
         integer lxss, lxyzs

         lxss  = 1
         lxyzs = 1

         sz = 3*lxss*2*lelv + lxyzs*lelv

         n = sz

         allocate(cb_fastg(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_fastg$',ierr)
         cb_fastg = 0

      end subroutine init
      end module fastg_mod
