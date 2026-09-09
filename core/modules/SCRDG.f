      module scrdg_mod

      implicit none

      real, allocatable, target :: cb_scrdg(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        conduct.f (fwght: uf(lf))
         integer sz_fwght
c        convect.f (conv_rhs_dg_aliased: uf,uxf,uyf,uzf,upwind_wgt(lf))
         integer sz_aliased
c        convect.f (conv_rhs_dg: + beta_c,jaco_c(lx1*lz1),
c        beta_f,jaco_f,ufine(lxd*lzd))
         integer sz_rhs_dg
c        convect.f (conv_bdry_dg_weak/conv_rhs_dg_weak: + us(lf))
         integer sz_dg_weak
         integer lf

         lf = lx1*lz1*2*ldim*lelt

         sz_fwght   = lf
         sz_aliased = 5*lf
         sz_rhs_dg  = 5*lf + 2*(lx1*lz1) + 3*(lxd*lzd)
         sz_dg_weak = 6*lf + 2*(lx1*lz1) + 3*(lxd*lzd)

         n = max(sz_fwght, sz_aliased, sz_rhs_dg, sz_dg_weak)

         allocate(cb_scrdg(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrdg$',ierr)
         cb_scrdg = 0

      end subroutine init
      end module scrdg_mod
