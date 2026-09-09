      module scrmg_mod

      implicit none

      real, allocatable, target :: cb_scrmg(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrmg
         integer sz_lelt4   ! hmholtz.f/hsmg.f/navier6.f (4 planes)
         integer sz_avg     ! prepost.f (tax: LX1*LY1*LELT*LDIMT)
         integer sz_geom    ! bdry.f nekasgn / subs1.f: xrm1..zsm1 (6 face planes)
         integer sz_edge    ! navier8.f: edge/enum/fnum tables (24+12+6)*lelt

         sz_lelt4 = 4*lx1*ly1*lz1*lelt
         sz_avg   = lx1*ly1*lelt*ldimt
         sz_geom  = 6*lx1*ly1
         sz_edge  = 42*lelt

         nscrmg = max(sz_lelt4, sz_avg, sz_geom, sz_edge)

         allocate(cb_scrmg(nscrmg), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrmg$',ierr)
         cb_scrmg = 0

      end subroutine init
      end module scrmg_mod
