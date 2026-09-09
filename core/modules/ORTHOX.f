      module orthox_mod

      implicit none

      real, allocatable, target :: cb_orthox(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        navier4.f (econj: pbar,pnew,pbrr(lx2*ly2*lz2*lelv)); dominates
c        every other ORTHOX site (setrhs/gensoln/updtset/updrhse,
c        induct.f's setrhsp/gensolnp/econjp: only pbar,pnew)
         integer sz

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

         sz = 3*lx2*ly2*lz2*lelv

         n = sz

         allocate(cb_orthox(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_orthox$',ierr)
         cb_orthox = 0

      end subroutine init
      end module orthox_mod
