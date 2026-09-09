      module cbplan_vol_ms_mod

      implicit none

      real, allocatable, target :: cb_cbplan_vol_ms(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        multimesh.f (plan3_vol_ms/plan4_vol_ms:
c        vxcp,dvxc,vycp,dvyc,vzcp,dvzc(lx1*ly1*lz1*lelv) x6,
c        resbc(lx1*ly1*lz1*lelv,ldim+1))
         integer sz

         logical icalld
         save    icalld
         data    icalld /.false./

         if (icalld) return
         icalld = .true.

         sz = (7+ldim)*lx1*ly1*lz1*lelv

         n = sz

         allocate(cb_cbplan_vol_ms(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_cbplan_vol_ms$',ierr)
         cb_cbplan_vol_ms = 0

      end subroutine init
      end module cbplan_vol_ms_mod
