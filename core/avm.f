c---------------------------------------------------------------
      subroutine setdefault_avm
      include 'SIZE'
      include 'AVM'

      do i=1,ldimt
        avm_ncut(i) = 1
        avm_c1(i) = 1.0
        ifcont(i) = .false.
        ifupwindavm(i) = .false.
      enddo

      return
      end
c---------------------------------------------------------------
      real function avm_vdiff(ix,iy,iz,e,cx,cy,cz)
c
c c1 and ncut a user tuneable control parameters. 
c Set c1 = 1.0 and reduce/increase as much possible/required,
c depending on your application.
c Typically ncut 1 or 2 works well. 
c Note, avoid using lx1 < 6! 
c
      include 'SIZE'
      include 'TOTAL'
      include 'AVM'
      include 'LVLSET'

      integer ix, iy, iz, e
      real c1, c2
      integer ncut

      real visc(lx1,ly1,lz1,lelt)
      save visc

      common /SCRMG/ r (lx1*ly1*lz1,lelt),
     $               tx(lx1*ly1*lz1,lelt),
     $               ty(lx1*ly1*lz1,lelt),
     $               tz(lx1*ly1*lz1,lelt)

      parameter (lm=40)
      parameter (lm2=lm*lm)
      real hpf_filter(lm2)
      real hpf_op(lx1*lx1,ldimt1)
      save hpf_op

      integer ibuild(ldimt1)
      save    ibuild

      save    icalld
      data    icalld / 0 /

      real h0,h0max
      real viscc(8,lelt)
      real cx(lx1,ly1,lz1,1) 
      real cy(lx1,ly1,lz1,1) 
      real cz(lx1,ly1,lz1,1)

      if (ix*iy*iz*e .ne. 1) then ! use cache
         avm_vdiff = max(1e-10,visc(ix,iy,iz,e))
         return
      endif

      if (icalld.eq.0) then
         iffilter(ifield) = .false.
         do i = 1,ldimt1
            ibuild(i) = 0
         enddo
         icalld = 1
      endif

      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv
      c2   = 0.5

      ! compute residual
      if (ifield.eq.1) then
         if (nid.eq.0) write(6,*) 'avm not supported for ifield=1 !'
         call exitt
      else
         if (ibuild(ifield).eq.0) then
           if(nid.eq.0)then
             write(6,*)'AVM Parameters for ifield:',ifield
             write(6,*)'ncut:',avm_ncut(ifield-1)
             write(6,*)'c1:',avm_c1(ifield-1)
             write(6,*)'C0 continuous:', ifcont(ifield-1)
             write(6,*)
           endif

           call hpf_trns_fcn(hpf_filter,avm_ncut(ifield-1))
           call build_hpf_mat(hpf_op(1,ifield),hpf_filter,.false.)
           ibuild(ifield) = ibuild(ifield) + 1
         endif
         call build_hpf_fld(r,t(1,1,1,1,ifield-1),hpf_op(1,ifield),
     $                      lx1,lz1)

         psave = param(99)
         param(99) = 0
         call copy(tx,r,n)
         if(ifield.eq.ifld_clsr)then
           call conv_clsr(r,tx)
         elseif(ifield.eq.ifld_tlsr)then
           call conv_tlsr(r,tx)
         else
           call convop(r,tx)
         endif
         param(99) = psave

         ! normalize
         uavg = gl2norm(t(1,1,1,1,ifield-1),n)
         call cfill(tx,uavg,n)
         call sub2(tx,t(1,1,1,1,ifield-1),n)
         uinf = 1.
         if(uavg.gt.0) uinf = glamax(tx, n)
      endif

      ! evaluate arificial viscosity
      h0max = 1e-10
      do ie=1,nelv
        h0max = max(h0max,deltaf(1,1,1,ie))
      enddo

      uinf = 1./uinf
      do ie = 1,nelv
      do i = 1,nxyz
         h0 = deltaf(i,1,1,ie)
         if(if3d)then
           vmax = sqrt(cx(i,1,1,ie)**2 + cy(i,1,1,ie)**2 
     $               + cz(i,1,1,ie)**2)
         else
           vmax = sqrt(cx(i,1,1,ie)**2 + cy(i,1,1,ie)**2)
         endif
         vismax = c2 * h0max * vmax
         visc(i,1,1,ie) = min(vismax, 
     $        avm_c1(ifield-1)*h0**2 * abs(r(i,ie))*uinf)
      enddo
      enddo

      ! make it piecewise constant
      do ie = 1,nelv
         vmax = vlmax(visc(1,1,1,ie),nxyz)
         call cfill(visc(1,1,1,ie),vmax,nxyz)
      enddo

      ! make it P1 continuous
      if (ifcont(ifield-1)) then
         call dsop (visc,'max',lx1,ly1,lz1)
         do ie = 1,nelv 
           viscc(1,ie) = visc(1  ,1  ,1  ,ie)
           viscc(2,ie) = visc(lx1,1  ,1  ,ie)
           viscc(3,ie) = visc(1  ,ly1,1  ,ie)
           viscc(4,ie) = visc(lx1,ly1,1  ,ie)

           viscc(5,ie) = visc(1  ,1  ,lz1,ie)
           viscc(6,ie) = visc(lx1,1  ,lz1,ie)
           viscc(7,ie) = visc(1  ,ly1,lz1,ie)
           viscc(8,ie) = visc(lx1,ly1,lz1,ie)
         enddo 
         call map_c_to_f_h1_bilin(visc,viscc)
      endif

      if (mod(istep,10).eq.0) then
        vismx = glamax(visc,n)
        vismn = glamin(visc,n)
        visav = glsc2(visc,bm1,n)/volvm1
        if (nio.eq.0) write(6,10) time,vismx,vismn,visav,ifield
 10     format(1p4e12.4,' AVM',i6)
      endif

      avm_vdiff = visc(ix,iy,iz,e)

      return
      end
c---------------------------------------------------------------
      real function deltaf(ix,iy,iz,iel)
c
c     compute characteristic sgs length scale
c     defined as min of GLL nodes within an elements but many
c     others possible.
c
      include 'SIZE'
      include 'TOTAL'

      real dx(lx1,ly1,lz1,lelt)
      save dx 

      data icalld /0/
      save icalld

      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv

      if (icalld.eq.0 .or. ifmvbd) then
         dinv = 1./ldim
         do ie = 1,nelv
c            volavg = 0
c            do i  = 1,nxyz
c               volavg = volavg + bm1(i,1,1,ie)
c            enddo
c            dd = (volavg**dinv)/(lx1-1)

            dd = dxmin_e(ie)

            call cfill(dx(1,1,1,ie),dd,nxyz) 
         enddo
         icalld = 1
      endif

      deltaf = dx(ix,iy,iz,iel)

      end 
c---------------------------------------------------------------
      subroutine axhelm_avm (au,u,imsh,isd)
      implicit none  
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'CTIMER'
      include 'TSTEP'
      include 'AVM'
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer imsh,isd,iz,nel,ntot,nxy,nyz,nxz,nxyz
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)
      real dudr,duds,dudt,tmp1,tmp2,tmp3

      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

      integer e

      common /avmtemp/ avmau(lx1,ly1,lz1,lelt)
      real avmau

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nel=nelt
      if (imsh.eq.1) nel=nelv

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NEL

      call rzero(avmau,ntot)

      do 100 e=1,nel
        IF (ldim.EQ.2) THEN
           call mxm  (dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           call mxm  (u(1,1,1,e),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
           endif
           call col2 (tmp1,avm_diff(1,1,1,e),nxyz)
           call col2 (tmp2,avm_diff(1,1,1,e),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (avmau(1,1,1,e),tm1,nxyz)
           call add2 (avmau(1,1,1,e),tm2,nxyz)
        else
           call mxm(dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,e),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,e),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
           call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
           endif
           call col2 (tmp1,avm_diff(1,1,1,e),nxyz)
           call col2 (tmp2,avm_diff(1,1,1,e),nxyz)
           call col2 (tmp3,avm_diff(1,1,1,e),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
           call add2 (avmau(1,1,1,e),tm1,nxyz)
           call add2 (avmau(1,1,1,e),tm2,nxyz)
           call add2 (avmau(1,1,1,e),tm3,nxyz)
        endif
C
 100  continue

      call add2(au,avmau,ntot)

      return
      end
