      module c_is1_mod

      implicit none

      integer*8, allocatable, target :: cb_c_is1(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr
         integer n
c        connect1.f (setup_topo/setup_mesh_dssum: glo_num(lx1*ly1*lz1*lelv))
         integer sz_connect1
c        drive1.f (nek_init: glo_num(lx1*ly1*lz1,lelt))
         integer sz_drive1
c        dssum.f (gtpp_gs_setup: glo_num(lx1,ly1,lz1,lelv))
         integer sz_dssum_v
c        dssum.f (gs_setup_ms: glo_num(lx1*ly1*lz1*lelt))
         integer sz_dssum_t
c        gmres.f (set_overlap2: glo_num(lxs*lys*lzs*lelv), lxs=lys=lzs=1)
         integer sz_gmres
c        convect.f (setup_dg_gs: glo_num_face(lf),
c        glo_num_vol((lx1+2)*(ly1+2)*(lz1+2)*lelt)); dominates
         integer sz_convect
c        hsmg.f (hsmg_setup_dssum: glo_num((lx1+2)*(ly1+2)*(lz1+2)*lelv))
         integer sz_hsmg_v
c        hsmg.f (h1mg_setup_dssum: glo_num((lx1+2)*(ly1+2)*(lz1+2)*lelt))
         integer sz_hsmg_t
c        hrefine.f (h_refine_usrdat2: glo_num(lx1*ly1*lz1*lelt))
         integer sz_hrefine
         integer lf, lxyzp

         lf    = lx1*lz1*2*ldim*lelt
         lxyzp = (lx1+2)*(ly1+2)*(lz1+2)

         sz_connect1 = lx1*ly1*lz1*lelv
         sz_drive1   = lx1*ly1*lz1*lelt
         sz_dssum_v  = lx1*ly1*lz1*lelv
         sz_dssum_t  = lx1*ly1*lz1*lelt
         sz_gmres    = lelv
         sz_convect  = lf + lxyzp*lelt
         sz_hsmg_v   = lxyzp*lelv
         sz_hsmg_t   = lxyzp*lelt
         sz_hrefine  = lx1*ly1*lz1*lelt

         n = max(sz_connect1, sz_drive1, sz_dssum_v, sz_dssum_t,
     $           sz_gmres, sz_convect, sz_hsmg_v, sz_hsmg_t, sz_hrefine)

         allocate(cb_c_is1(n), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_c_is1$',ierr)
         cb_c_is1 = 0_8

      end subroutine init
      end module c_is1_mod
