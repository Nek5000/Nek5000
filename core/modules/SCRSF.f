      module scrsf_mod

      implicit none

      real, allocatable, target :: cb_scrsf(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrsf
         integer sz_stress   ! subs1.f (stress: t11,t22,t33,hii)
         integer sz_avg      ! navier5.f (auto_averager: ta(lt,ldimt))

         sz_stress = 4*lx1*ly1*lz1*lelt
         sz_avg    = lx1*ly1*lz1*lelt*ldimt

         nscrsf = max(sz_stress, sz_avg)

         allocate(cb_scrsf(nscrsf), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrsf$',ierr)
         cb_scrsf = 0

      end subroutine init
      end module scrsf_mod
