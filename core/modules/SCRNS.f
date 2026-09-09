      module scrns_mod

      implicit none

      real, allocatable, target :: cb_scrns(:)
      real, allocatable, target :: cb_resdmp(:)  ! restart T/ps dump (was common /cbresdmp/); separate from cb_scrns

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrns

         integer sz_induct
         integer sz_postpro
         integer sz_lelt7
         integer sz_lelv9
         integer sz_navier8
         integer sz_convect
         integer sz_prepost
         integer sz_hrefine

         sz_induct  = 14*lxd*lyd*lzd + 6*lx1*ly1*lz1
         sz_postpro = 60000
         sz_lelt7   = 7*lx1*ly1*lz1*lelt
         sz_lelv9   = 9*lx1*ly1*lz1*lelv + lx2*ly2*lz2*lelv
         sz_navier8 = 432*lelt
         sz_convect = lxd*lyd*lzd*lelv*ldim
         sz_prepost = 1 + 3*lxo*lxo*lxo*lelt
         sz_hrefine = max(512,lelt)*(3*lx1*ly1*lz1 + 2*lx1*lx1)

         nscrns = max(sz_induct, sz_postpro, sz_lelt7,
     $                sz_lelv9, sz_navier8, sz_convect, sz_prepost,
     $                sz_hrefine)

         allocate(cb_scrns(nscrns), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrns$',ierr)
         cb_scrns = 0

         allocate(cb_resdmp(lx1*ly1*lz1*lelt*ldimt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_resdmp$',ierr)
         cb_resdmp = 0

      end subroutine init
      end module scrns_mod
