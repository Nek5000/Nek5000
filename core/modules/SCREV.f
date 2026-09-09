      module screv_mod

      implicit none

      real, allocatable, target :: cb_screv(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrev

c        SCREV is a shared scratch pool carved into different partitions by its
c        consumers; size it for the LARGEST one. Consumers:
c          two-plane family : 2*lx1*ly1*lz1*lelt (dominant in 3D)
c            (drive2.f SII/SIII, eigsolv.f X1/Y1, subs1.f dpc/p1,
c             subs2.f ei2/ei3, navier1.f H1/H2 & RKZ1/RKZ2, navier6.f x(2*ltotd))
c          subs2.f STSMASK/UPDMSYS hfmask+hvmask : 6*lx1*lz1*lelt + lx1*ly1*lz1*lelt
c          postpro.f filter_s0 wk1+wk2           : lx1*lx1*lx1*lelt + lx1*lx1*lx1
         nscrev = max( 2*lx1*ly1*lz1*lelt
     $               , 6*lx1*lz1*lelt + lx1*ly1*lz1*lelt
     $               , lx1*lx1*lx1*lelt + lx1*lx1*lx1 )

         allocate(cb_screv(nscrev), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_screv$',ierr)
         cb_screv = 0

      end subroutine init
      end module screv_mod
