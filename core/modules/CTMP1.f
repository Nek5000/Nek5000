      module ctmp1_mod

      implicit none

      real, allocatable, target :: cb_ctmp1(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nctmp1

c        LX1*LY1*LZ1*LELT family max (coef.f, subs1.f, subs2.f: 4
c        "planes"); covers every other same-scaled site (<=4 planes)
         integer sz_lelt4
c        LX1*LY1*LZ1*LELV family max (navier1.f RKX1..RKX4: 4
c        "planes"); covers every other same-scaled site (<=4 planes)
         integer sz_lelv4
c        navier1.f (local_grad3_t/local_grad2_t: ur..wt, 9 fields);
c        single-element scratch, no LELT/LELV factor
         integer sz_navier1
c        convect.f (dealiasing scratch, ur,us,ut,ju,ud,tu); LXD can
c        exceed LX1 by up to 1.5x, so this can exceed sz_navier1
         integer sz_convect
c        genxyz.f (H,XCRVED,YCRVED,ZCRVED,ZGML,WORK); single-element
         integer sz_genxyz
c        hmholtz.f (DUDR,DUDS,DUDT,TMP1,TMP2,TMP3); single-element
         integer sz_hmholtz
c        navier5.f (filter_d2: w(lt,lelt)+ur+us+ut+w1(2lt))
         integer sz_navier5
c        navier6.f (set_overlap: mask, integer, packed storage)
         integer sz_navier6
c        ic.f (restart: TDUMP(LXYZR,LPSC9), real*4, LXR=LX1+6 padding
c        for restarting from a different mesh order)
         integer sz_icrestart
c        prepost.f (check_ioinfo/outfld: TDUMP(LXYZ,LPSC9), real*4)
         integer sz_prepost
c        reader_re2.f (readp_re2_curve: vi(li,nrmax), integer)
         integer sz_re2vi
c        bdry.f (TRST3D: DRM1,DRTM1,DSM1,DSTM1,WGS)
         integer sz_bdry

         sz_lelt4     = 4*lx1*ly1*lz1*lelt
         sz_lelv4     = 4*lx1*ly1*lz1*lelv
         sz_navier1   = 9*lx1*ly1*lz1
         sz_convect   = 6*lxd*lyd*lzd
         sz_genxyz    = 12*lx1 + 3*lx1*lx1
         sz_hmholtz   = 6*lx1*ly1*lz1
         sz_navier5   = lx1*ly1*lz1*lelt + 5*lx1*ly1*lz1
         sz_navier6   = 2*lx1*ly1*lz1*lelt
         sz_icrestart = ((lx1+6)*(ly1+6)*(lz1+6)*(ldimt+9) + 1)/2
         sz_prepost   = (lx1*ly1*lz1*(ldimt+10) + 1)/2
         sz_re2vi     = 102*lelt
         sz_bdry      = 2*lx1*lx1 + 3*lx1*ly1

         nctmp1 = max(sz_lelt4, sz_lelv4, sz_navier1, sz_convect,
     $                sz_genxyz, sz_hmholtz, sz_navier5, sz_navier6,
     $                sz_icrestart, sz_prepost, sz_re2vi, sz_bdry)

         allocate(cb_ctmp1(nctmp1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_ctmp1$',ierr)
         cb_ctmp1 = 0

      end subroutine init
      end module ctmp1_mod
