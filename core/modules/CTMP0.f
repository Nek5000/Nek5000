      module ctmp0_mod

      implicit none

      real, allocatable, target :: cb_ctmp0(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nctmp0
         integer lxm1

c        LX1*LY1*LZ1*LELT family max (coef.f, hmholtz.f, navier4.f,
c        subs1.f, subs2.f, navier6.f: 2 "planes")
         integer sz_lelt2
c        LX1*LY1*LZ1*LELV family max (subs1.f: 2 "planes")
         integer sz_lelv2
c        LX2*LY2*LZ2*LELV family max (navier1.f, perturb.f: 2 "planes")
         integer sz_lelv2x2
c        navier1.f (op_curl-style dealiasing scratch: TA1D,DFID);
c        LXD can exceed LX1 by up to 1.5x
         integer sz_dealias
c        hmholtz.f (h1mg fdm setup: W,aa,bb,mask)
         integer sz_hmholtz_a
c        hmholtz.f (hxdg scratch: w)
         integer sz_hmholtz_b
c        convect.f (grad_rstd/map_faced dealiasing scratch: w)
         integer sz_convect
c        hsmg.f (hsmg_setup_dbg: w(100000), fixed); also safely
c        covers every small fixed-size CTMP0 site in the source
         integer sz_hsmg_big
c        ic.f (mapab-style restart interp: xa,xb,xc,zgmr,wgtr; LXR=
c        LX1+6 padding for restarting from a different mesh order)
         integer sz_icrestart
c        navier5.f (q_filter: intw,intt,wk1,wk2,zgmv,wgtv,zgmp,wgtp,
c        tmax,omax)
         integer sz_navier5_a
c        navier5.f (comp_vort3: ur..wt, 9 fields, no LELT/LELV factor)
         integer sz_navier5_b
c        connect2.f (inrtia: VO,XYZI,CG,TI,WORK; LELT-linear, no LX1)
         integer sz_connect2
c        gmres.f (wk1,wk2)
         integer sz_gmres
c        navier8.f (map_c_to_f_h1_bilin/map_f_to_c_h1_bilin: w,v)
         integer sz_navier8
c        hsmg.f (hsmg_setup_fast1d: b,w; LXM=LX1+2)
         integer sz_hsmg_lxm
c        map2.f (find_con: eid8(4*lelt),vtx8(lelt*(2**ldim+1)))
         integer sz_map2
c        genxyz.f (xyzquad: w(4*lx1*ly1*lz1,2),zg(3))
         integer sz_genxyz

         lxm1 = lx1 + 2

         sz_lelt2      = 2*lx1*ly1*lz1*lelt
         sz_lelv2      = 2*lx1*ly1*lz1*lelv
         sz_lelv2x2    = 2*lx2*ly2*lz2*lelv
         sz_dealias    = lx1*ly1*lz1*lelv + 2*lxd*lyd*lzd*lelv
         sz_hmholtz_a  = 3*lx1*lx1 + lx1*ly1*lz1*lelt
         sz_hmholtz_b  = 4*lx1*lz1*ldim*lelt
         sz_convect    = 2*(2*lxd)**ldim
         sz_hsmg_big   = 100000
         sz_icrestart  = 2*(lx1+6)*(ly1+6)*(lz1+6)
     $                 + lx1*ly1*(lz1+6) + 2*(lx1+6)
         sz_navier5_a  = 2*lx1*lx1 + lx1*lx1*lx1*lelt + lx1*lx1*lx1
     $                 + 4*lx1 + 203
         sz_navier5_b  = 9*lx1*ly1*lz1
         sz_connect2   = 7*lelt + 12
         sz_gmres      = 2*lgmres
         sz_navier8    = 2*lx1*lx1 + 2*lx1*(ldim-1)*lelt
         sz_hsmg_lxm   = 4*lxm1*lxm1
         sz_map2       = 5*lelt + lelt*(2**ldim)
         sz_genxyz     = 8*lx1*ly1*lz1 + 3

         nctmp0 = max(sz_lelt2, sz_lelv2, sz_lelv2x2, sz_dealias,
     $                sz_hmholtz_a, sz_hmholtz_b, sz_convect,
     $                sz_hsmg_big, sz_icrestart, sz_navier5_a,
     $                sz_navier5_b, sz_connect2, sz_gmres,
     $                sz_navier8, sz_hsmg_lxm, sz_map2, sz_genxyz)

         allocate(cb_ctmp0(nctmp0), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_ctmp0$',ierr)
         cb_ctmp0 = 0

      end subroutine init
      end module ctmp0_mod
