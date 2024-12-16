c---------------------------------------------------------------------
      subroutine setdefault_svv
      implicit none
      include 'SIZE'
      include 'SVV'

      integer i

      svv_c0_fluid = 0.1
      svvcut_fluid = (lx1-1.0)/2.0

      do i=1,ldimt+1
        svv_c0(i) = 0.1
        svvcut(i) = (lx1-1.0)/2.0
        ifnlsvv(i) = .false.
        ifupwindsvv(i) = .false.

        svv_k0(i) = 1.5
        svv_k1(i) = 1.5
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine setmu_svv(phi,cx,cy,cz)
c      
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      real phi(lx1,ly1,lz1,lelt)
      real csf(lx1,ly1,lz1,lelt)
c
      integer icalld
      save icalld
      data icalld /0/

      real cx(1),cy(1),cz(1)
c
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv
      
      if(icalld .eq. 0)then
        call getconvD
         call copy(svvf,jacm1,ntot)
         do i=1,ntot
            if(if3d)then
               svvf(i,1,1,1) = svvf(i,1,1,1)**(1.0/3.0)
            else
               svvf(i,1,1,1) = svvf(i,1,1,1)**(1.0/2.0)
            endif
         enddo
         call cmult(svvf,2.0,ntot)
         
c     Scale with mu_0/N
         call cmult(svvf,1.0/(lx1-1.0),ntot)
         icalld = 1
      endif

c     Scale with advection velocity
      call col3(svvmu,cx,cx,ntot)
      call addcol3(svvmu,cy,cy,ntot)
      if(if3d) call addcol3(svvmu,cz,cz,ntot)
      do i=1,ntot
         svvmu(i,1,1,1) = svvmu(i,1,1,1)**0.5
      enddo
      
      call col2(svvmu,svvf,ntot)
      call cmult(svvmu,svv_c0(ifield),ntot)
      call copy(svvprec,svvmu,ntot)
      if(ifnlsvv(ifield))then
        call getnlsvvsf(phi,csf)
        call col2(svvmu,csf,ntot)
      endif
      
      return
      end
c---------------------------------------------------------------------
      subroutine getnlsvvsf(phi,csf)
c      
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      real phi(lx1,ly1,lz1,lelt)
      real csf(lx1,ly1,lz1,lelt)
      real cdi(lx1,ly1,lz1,lelt)
      
      real gradx(lx1,ly1,lz1,lelt)
      real grady(lx1,ly1,lz1,lelt)
      real gradz(lx1,ly1,lz1,lelt)

      real switch(lelt)
c
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      call gradm1(gradx,grady,gradz,phi)
      call opcolv(gradx,grady,gradz,bm1)
      call opdssum(gradx,grady,gradz)
      call opcolv(gradx,grady,gradz,binvm1)
      
      if(if3d)then
         call vdot3(cdi,gradx,grady,gradz,gradx,grady,gradz,ntot)
      else
         call vdot2(cdi,gradx,grady,gradx,grady,ntot)
      endif
      
c      call col2(cdi,bm1,ntot)
      
      call cmult(cdi,1.0/glamax(cdi,ntot),ntot)

      call get_visc_switch(phi,switch)
            
      do ie=1,nelv
         do i=1,nxyz
           cnl = (lx1-1.0)**svv_k1(ifield)
            csf(i,1,1,ie) = switch(ie)*max(0.5,cnl*cdi(i,1,1,ie))            
         enddo
      enddo
c     
      return
      end
c---------------------------------------------------------------------
      subroutine get_visc_switch(phi,switch)
c      
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      real phi(lx1,ly1,lz1,lelt)
      real switch(lelt)
      real eratio
c      
      real phim(lx1,ly1,lz1,lelt)
      real tmph(lx1,ly1,lz1)
      real tmpl(lx1,ly1,lz1)
      real suml
c      
      real maskh(lx1,ly1,lz1)
      save maskh
      
      real maskl(lx1,ly1,lz1)
      save maskl

      real threshold
      save threshold
c
      integer icalld
      save icalld
      data icalld /0/
c
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv
      
      call getspec3D(phi,phim)
      
      if(icalld .eq. 0)then
         call rzero(maskh,nxyz)
         call rzero(maskl,nxyz)
         
         if(if3d)then
            do k=1,lz1
               do j=1,ly1
                  do i=1,lx1
                     if(i.eq.1 .XOR. j.eq.1 .XOR. k.eq.1)maskl(i,j,k)=1.
                     if(i.eq.lx1 .XOR. j.eq.ly1 .XOR. k.eq.lz1)
     $                    maskh(i,j,k)=1.
                  enddo
               enddo
            enddo
         else
            do j=1,ly1
               do i=1,lx1
                  if(i.eq.1 .XOR. j.eq.1)maskl(i,j,1)=1.
                  if(i.eq.lx1 .XOR. j.eq.ly1) maskh(i,j,1)=1.
               enddo
            enddo
         endif
         threshold = 10.0**((lx1-1.)/svv_k0(ifield))      
         icalld = 1
      endif
     
      do ie=1,nelv
         call col3(tmph,maskh,phim(1,1,1,ie),nxyz)
         call col3(tmpl,maskl,phim(1,1,1,ie),nxyz)

         suml = vlsum(tmpl,nxyz)

         if(suml .gt. 0.)then
            eratio = vlsum(tmph,nxyz)/suml
         else
            eratio = vlsum(tmph,nxyz)
         endif

         switch(ie) = visc_switch(eratio*threshold)
      enddo
     
      return
      end
c---------------------------------------------------------------------
      real function visc_switch(a)
      real a
      real PI

      PI = 2.0*acos(0.0)
      
      if(a.ge.1.)then
         visc_switch = 1.
      else
         visc_switch = abs(0.5*(2.0*a + (1./PI)*sin(PI*(2.0*a-1.0))))
         visc_switch = max(1.0e-15,visc_switch)
      endif

      return
      end
c---------------------------------------------------------------------      
      subroutine getspec3D(phi,phim)
c     
      include 'SIZE'
      include 'TOTAL'
c
      real phi(lx1,ly1,lz1,lelt)
      real phim(lx1,ly1,lz1,lelt)
c      
      real bl(lx1,lx1)
      real blt(lx1,lx1)
c      
      real bltinv(lx1,lx1)
      save bltinv

      real blinv(lx1,lx1)
      save blinv

      real tmp1(lx1*ly1*lz1)
      real tmp2(lx1*ly1*lz1)
c
      integer icalld
      save icalld
      data icalld /0/
c
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv
      
      if(icalld .eq. 0)then
         call leg1D(bl,blinv,blt,bltinv)
         icalld = 1
      endif

      do ie=1,nelv
         if(.not.if3d)then
            call mxm(blinv,lx1,phi(1,1,1,ie),lx1,tmp1,lx1)
            call mxm(tmp1,lx1,bltinv,lx1,phim(1,1,1,ie),lx1)
         else
            call mxm(blinv,lx1,phi(1,1,1,ie),lx1,tmp1,lx1*lx1)
            ix = 1
            iy = 1
            do iz=1,lz1
               call mxm(tmp1(ix),lx1,bltinv,lx1,tmp2(iy),lx1)
               ix = ix+lx1*lx1
               iy = iy+lx1*lx1
            enddo
            call mxm(tmp2,lx1*lx1,bltinv,lx1,phim(1,1,1,ie),lx1)
         endif
      enddo

      do i=1,ntot
         phim(i,1,1,1) = abs(phim(i,1,1,1))
      enddo
      return
      end      
c---------------------------------------------------------------------
      subroutine getconvD
c
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      real bl(lx1,lx1), blinv(lx1,lx1)
      real blt(lx1,lx1), bltinv(lx1,lx1)
      real qdiag(lx1,lx1)
      real tmp(lx1,lx1)
      real cmat(lx1,lx1), cmatT(lx1,lx1)
c
      call leg1D(bl,blinv,blt,bltinv)
      call diffFilter1D(qdiag)

      call mxm(bl,lx1,qdiag,lx1,tmp,lx1)
      call mxm(tmp,lx1,blinv,lx1,cmat,lx1)

      call mxm(bltinv,lx1,qdiag,lx1,tmp,lx1)
      call mxm(tmp,lx1,blt,lx1,cmatT,lx1)

c     Convoluted derivatives
      call mxm(cmat,lx1,dxm1,lx1,cdxm1,lx1)
      call mxm(dxtm1,lx1,cmatT,lx1,cdxtm1,lx1)

      call mxm(cmat,ly1,dym1,ly1,cdym1,ly1)
      call mxm(dytm1,ly1,cmatT,ly1,cdytm1,ly1)

      if(if3d)then
         call mxm(cmat,lz1,dzm1,lz1,cdzm1,lz1)
         call mxm(dztm1,lz1,cmatT,lz1,cdztm1,lz1)
      endif

c     Store original derivatives
      call copy(odxm1,dxm1,lx1*lx1)
      call copy(odxtm1,dxtm1,lx1*lx1)

      call copy(odym1,dym1,ly1*ly1)
      call copy(odytm1,dytm1,ly1*ly1)

      if(if3d)then
         call copy(odzm1,dzm1,lz1*lz1)
         call copy(odztm1,dztm1,lz1*lz1)
      endif
      
      return
      end
c---------------------------------------------------------------------      
      subroutine leg1D(bl,blinv,blt,bltinv)
c      
c     1D Vandermonde matrix corresponding to Legendre Basis
c      
      include 'SIZE'
      include 'TOTAL'
c
      real bl(lx1,lx1), blinv(lx1,lx1)
      real blt(lx1,lx1), bltinv(lx1,lx1)
      real rmult(lx1)
      integer indr(lx1),indc(lx1),ipiv(lx1)
     
      call LGLL(blinv,zgm1,lx1)
      call transpose(bl,lx1,blinv,lx1)

      call copy(blinv,bl,lx1*lx1)

      call gaujordf(blinv,lx1,lx1,indr,indc,ipiv,ierr,rmult)

      call transpose(blt,lx1,bl,lx1)
      call transpose(bltinv,lx1,blinv,lx1)
      
      return
      end
c---------------------------------------------------------------------
      subroutine LGLL (B,Z,NZ)
c     Get all the 1d Legendre basis at Z for polynomial order NZ-1
      real B(NZ,NZ)
      real Z(NZ)
      
      N = NZ-1

      if(NZ .eq. 1)then
         B(1,1) = 0.
         return
      endif

      do i=1,NZ
         do j=1,NZ
            B(i,j) = PNLEG(Z(j),i-1)
         enddo
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine diffFilter1D(q)
c      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      real q(lx1,lx1)
      real a,k
      real n
      integer i
c
      n = lx1-1.0

      call ident(q,lx1)
      
      do i=1,lx1
         k = i-1
         q(i,i) = (k/n)**(svvcut(ifield)/2.0)
      enddo
      
      if(nid .eq. 0)then
         write(*,*)"---------Modal parameters for SVV filter---------"
         write(*,*)"Field is",ifield
         write(*,*)(q(i,i),i=1,lx1)
         write(*,*)"-------------------------------------------------"
      endif

      
      return
      end
c---------------------------------------------------------------------
      subroutine axhelm_svv (au,u,imsh,isd)
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'CTIMER'
      include 'TSTEP'
      include 'SVV'
      include 'SOLN'
C
      COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)
      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)

      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

      integer e

      common /svvtemp/ svvau(lx1,ly1,lz1,lelt)
      real svvau

      common /svvtemp2/ gux(lx1,ly1,lz1,lelt),
     $                  guy(lx1,ly1,lz1,lelt),
     $                  guz(lx1,ly1,lz1,lelt),
     $                  gdot(lx1,ly1,lz1,lelt)

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nel=nelt
      if (imsh.eq.1) nel=nelv

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NEL

      CALL RZERO (SVVAU,NTOT)

      if(ifupwindsvv(ifield))then
        call gradsvv(gux,guy,guz,u)
        call svvbdryfix

        if(if3d)then
          call vdot3(gdot,gux,guy,guz,svvnx,svvny,svvnz,ntot)
        else
          call vdot2(gdot,gux,guy,svvnx,svvny,ntot)
        endif
      endif

      do 100 e=1,nel
        IF (ldim.EQ.2) THEN
          if(ifupwindsvv(ifield))then
            svmin = vlmin(svvmask(1,1,1,e),nxyz)
            if(svmin.eq.1.0)then
              call col3(dudr,svvnr(1,1,1,e),gdot(1,1,1,e),nxyz)
              call col3(duds,svvns(1,1,1,e),gdot(1,1,1,e),nxyz)
            else
              goto 1990
            endif
            goto 1991
          endif
1990      call mxm  (cdxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
          call mxm  (u(1,1,1,e),lx1,cdytm1,ly1,duds,ly1)

1991       call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
           endif
           call mxm  (cdxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,cdym1,ly1,tm2,ly1)
           call add2 (svvau(1,1,1,e),tm1,nxyz)
           call add2 (svvau(1,1,1,e),tm2,nxyz)
        else
          if(ifupwindsvv(ifield))then
            svmin = vlmin(svvmask(1,1,1,e),nxyz)
            if(svmin.eq.1.0)then
              call col3(dudr,svvnr(1,1,1,e),gdot(1,1,1,e),nxyz)
              call col3(duds,svvns(1,1,1,e),gdot(1,1,1,e),nxyz)
              call col3(dudt,svvnt(1,1,1,e),gdot(1,1,1,e),nxyz)
            else
              goto 1992
            endif
            goto 1993
          endif
1992       call mxm(cdxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,e),lx1,cdytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,e),nxy,cdztm1,lz1,dudt,lz1)

1993       call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
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
           call mxm  (cdxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,cdym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,cdzm1,lz1,tm3,lz1)
           call add2 (svvau(1,1,1,e),tm1,nxyz)
           call add2 (svvau(1,1,1,e),tm2,nxyz)
           call add2 (svvau(1,1,1,e),tm3,nxyz)
        endif
 100  continue

      call col2(svvau,svvmu,ntot)

      call add2(au,svvau,ntot)
      taxhm=taxhm+(dnekclock()-etime1)
      return
      end
c---------------------------------------------------------------------
      subroutine setprec_svv(dpcm1,imsh,isd)
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'SVV'

      REAL            DPCM1 (LX1,LY1,LZ1,1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      REAL            HELM1(lx1,ly1,lz1,1), HELM2(lx1,ly1,lz1,1)
      REAL YSM1(LY1)

      common /svvdp/ sdpcm1(lx1,ly1,lz1,lelt)
      real sdpcm1

      nel=nelt
      if (imsh.eq.1) nel=nelv

      ntot = nel*lx1*ly1*lz1
      
      CALL RZERO(SDPCM1,NTOT)
      DO 1000 IE=1,NEL
        DO 320 IQ=1,lx1
        DO 320 IZ=1,lz1
        DO 320 IY=1,ly1
        DO 320 IX=1,lx1
           SDPCM1(IX,IY,IZ,IE) = SDPCM1(IX,IY,IZ,IE) + 
     $                          G1M1(IQ,IY,IZ,IE) * CDXTM1(IX,IQ)**2
  320   CONTINUE
        DO 340 IQ=1,ly1
        DO 340 IZ=1,lz1
        DO 340 IY=1,ly1
        DO 340 IX=1,lx1
           SDPCM1(IX,IY,IZ,IE) = SDPCM1(IX,IY,IZ,IE) + 
     $                          G2M1(IX,IQ,IZ,IE) * CDYTM1(IY,IQ)**2
  340   CONTINUE
        IF (LDIM.EQ.3) THEN
           DO 360 IQ=1,lz1
           DO 360 IZ=1,lz1
           DO 360 IY=1,ly1
           DO 360 IX=1,lx1
              SDPCM1(IX,IY,IZ,IE) = SDPCM1(IX,IY,IZ,IE) + 
     $                             G3M1(IX,IY,IQ,IE) * CDZTM1(IZ,IQ)**2
  360      CONTINUE
C
C          Add cross terms if element is deformed.
C
           IF (IFDFRM(IE)) THEN
              DO 600 IY=1,ly1,ly1-1
              DO 600 IZ=1,lz1,max(1,lz1-1)
              SDPCM1(1,IY,IZ,IE) = SDPCM1(1,IY,IZ,IE)
     $            + G4M1(1,IY,IZ,IE) * CDXTM1(1,1)*CDYTM1(IY,IY)
     $            + G5M1(1,IY,IZ,IE) * CDXTM1(1,1)*CDZTM1(IZ,IZ)
              SDPCM1(lx1,IY,IZ,IE) = SDPCM1(lx1,IY,IZ,IE)
     $            + G4M1(lx1,IY,IZ,IE) * CDXTM1(lx1,lx1)*CDYTM1(IY,IY)
     $            + G5M1(lx1,IY,IZ,IE) * CDXTM1(lx1,lx1)*CDZTM1(IZ,IZ)
  600         CONTINUE
              DO 700 IX=1,lx1,lx1-1
              DO 700 IZ=1,lz1,max(1,lz1-1)
                 SDPCM1(IX,1,IZ,IE) = SDPCM1(IX,1,IZ,IE)
     $            + G4M1(IX,1,IZ,IE) * CDYTM1(1,1)*CDXTM1(IX,IX)
     $            + G6M1(IX,1,IZ,IE) * CDYTM1(1,1)*CDZTM1(IZ,IZ)
                 SDPCM1(IX,ly1,IZ,IE) = SDPCM1(IX,ly1,IZ,IE)
     $            + G4M1(IX,ly1,IZ,IE) * CDYTM1(ly1,ly1)*CDXTM1(IX,IX)
     $            + G6M1(IX,ly1,IZ,IE) * CDYTM1(ly1,ly1)*CDZTM1(IZ,IZ)
  700         CONTINUE
              DO 800 IX=1,lx1,lx1-1
              DO 800 IY=1,ly1,ly1-1
                 SDPCM1(IX,IY,1,IE) = SDPCM1(IX,IY,1,IE)
     $                + G5M1(IX,IY,1,IE) * CDZTM1(1,1)*CDXTM1(IX,IX)
     $                + G6M1(IX,IY,1,IE) * CDZTM1(1,1)*CDYTM1(IY,IY)
                 SDPCM1(IX,IY,lz1,IE) = SDPCM1(IX,IY,lz1,IE)
     $              + G5M1(IX,IY,lz1,IE)*CDZTM1(lz1,lz1)*CDXTM1(IX,IX)
     $              + G6M1(IX,IY,lz1,IE)*CDZTM1(lz1,lz1)*CDYTM1(IY,IY)
  800         CONTINUE
           ENDIF

        ELSE  ! 2D

           IZ=1
           IF (IFDFRM(IE)) THEN
              DO 602 IY=1,ly1,ly1-1
                 SDPCM1(1,IY,IZ,IE) = SDPCM1(1,IY,IZ,IE)
     $                + G4M1(1,IY,IZ,IE) * CDXTM1(1,1)*CDYTM1(IY,IY)
                 SDPCM1(lx1,IY,IZ,IE) = SDPCM1(lx1,IY,IZ,IE)
     $          + G4M1(lx1,IY,IZ,IE) * CDXTM1(lx1,lx1)*CDYTM1(IY,IY)
  602         CONTINUE
              DO 702 IX=1,lx1,lx1-1
                 SDPCM1(IX,1,IZ,IE) = SDPCM1(IX,1,IZ,IE)
     $                + G4M1(IX,1,IZ,IE) * CDYTM1(1,1)*CDXTM1(IX,IX)
                 SDPCM1(IX,ly1,IZ,IE) = SDPCM1(IX,ly1,IZ,IE)
     $          + G4M1(IX,ly1,IZ,IE) * CDYTM1(ly1,ly1)*CDXTM1(IX,IX)
  702         CONTINUE
           ENDIF

        ENDIF
 1000 CONTINUE

      call col2(SDPCM1,svvprec,ntot)
      call dssum(SDPCM1,lx1,ly1,lz1)
      call add2(DPCM1,SDPCM1,ntot)
C
      return
      END
c---------------------------------------------------------------------
      subroutine gradsvv(ux,uy,uz,u)
c
c     Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)
c
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SVV'
c
      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz,1),uy(lxyz,1),uz(lxyz,1),u(lxyz,1)

      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)

      integer e

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelt

      N = lx1-1
      do e=1,nelt
         if (if3d) then
            call local_grad3(ur,us,ut,u,N,e,cdxm1,cdxtm1)
            do i=1,lxyz
               ux(i,e) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                             + us(i)*sxm1(i,1,1,e)
     $                             + ut(i)*txm1(i,1,1,e) )
               uy(i,e) = jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                             + us(i)*sym1(i,1,1,e)
     $                             + ut(i)*tym1(i,1,1,e) )
               uz(i,e) = jacmi(i,e)*(ur(i)*rzm1(i,1,1,e)
     $                             + us(i)*szm1(i,1,1,e)
     $                             + ut(i)*tzm1(i,1,1,e) )
            enddo
         else
            if (ifaxis) call setaxdy (ifrzer(e))
            call local_grad2(ur,us,u,N,e,cdxm1,cdytm1)
            do i=1,lxyz
               ux(i,e) =jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                            + us(i)*sxm1(i,1,1,e) )
               uy(i,e) =jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                            + us(i)*sym1(i,1,1,e) )
            enddo
         endif
      enddo
c
      return
      end
C----------------------------------------------------------------------     
      subroutine svvbdryfix
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'

      integer ie,ifc,i
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz
      integer ntot,nxyz
      integer icalld
      save icalld
      data icalld /0/

      real w1,w2,w3,w4,w5
      common /SCRNS/
     & w1(lx1*ly1*lz1*lelv)
     &,w2(lx1*ly1*lz1*lelv)
     &,w3(lx1*ly1*lz1*lelv)
     &,w4(lx1*ly1*lz1*lelv)
     &,w5(lx1*ly1*lz1*lelv)

      real ywd(lx1,ly1,lz1,lelt)
      real dxmax_e
      real dmax

      if(icalld.eq.0)then
        call distf(ywd,1,'O  ',w1,w2,w3,w4,w5)

        ntot = lx1*ly1*lz1*nelt
        nxyz = lx1*ly1*lz1

        dmax = dxmax_e(1)
        do i=1,ntot
          svvmask(i,1,1,1) = tanh(ywd(i,1,1,1)/(dmax))
        enddo

        icalld = 1
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine setUpwindSVVAVM(cx,cy,cz)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'

      real cx(1),cy(1),cz(1)

      common /svvtemp3/ cmag(lx1*ly1*lz1*lelt)
      real cmag

      integer ntot,i

      ntot = lx1*ly1*lz1*nelt

      call col3(cmag,cx,cx,ntot)
      call addcol3(cmag,cy,cy,ntot)
      if(if3d) call addcol3(cmag,cz,cz,ntot)
      call vsqrt(cmag,ntot)

      do i=1,ntot
        if(cmag(i).gt.0.0)then
          svvnx(i,1,1,1) = cx(i)/cmag(i)
          svvny(i,1,1,1) = cy(i)/cmag(i)
          if(if3d)svvnz(i,1,1,1) = cz(i)/cmag(i)
        else
          svvnx(i,1,1,1) = 0.0
          svvny(i,1,1,1) = 0.0
          svvnz(i,1,1,1) = 0.0
        endif
      enddo

      call vector_to_rst(svvnx,svvny,svvnz,svvnr,svvns,svvnt)

      return
      end
c---------------------------------------------------------------------
      subroutine axhelm_svv_fluid (au,u,imsh,isd)
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'CTIMER'
      include 'TSTEP'
      include 'SVV'
      include 'SOLN'
C
      COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)
      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)

      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

      integer e

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nel=nelt
      if (imsh.eq.1) nel=nelv

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NEL

      CALL RZERO (AU,NTOT)

      do 100 e=1,nel
        IF (ldim.EQ.2) THEN
          call mxm  (cdxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
          call mxm  (u(1,1,1,e),lx1,cdytm1,ly1,duds,ly1)

          call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
          call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
           endif
           call mxm  (cdxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,cdym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)
        else
          call mxm(cdxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,e),lx1,cdytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,e),nxy,cdztm1,lz1,dudt,lz1)
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
           call mxm  (cdxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,cdym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,cdzm1,lz1,tm3,lz1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)
           call add2 (au(1,1,1,e),tm3,nxyz)
        endif
 100  continue

      call col2(au,svvmu,ntot)

      call col2(au,bm1,ntot)
      call dssum(au,lx1,ly1,lz1)
      call col2(au,binvm1,ntot)

      taxhm=taxhm+(dnekclock()-etime1)
      return
      end
c---------------------------------------------------------------------
      subroutine svv_fluid(ix,iy,iz,e,svvx,svvy,svvz)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /svvfluid/ auvx(lx1,ly1,lz1,lelv),
     $                  auvy(lx1,ly1,lz1,lelv),
     $                  auvz(lx1,ly1,lz1,lelv),
     $                  dummy(lx1,ly1,lz1,lelv)
      real auvx,auvy,auvz

      real svvx,svvy,svvz
      integer ix,iy,iz,e

      integer ntot
      real dummy

      if(ix*iy*iz*e .eq. 1)then
        ntot = lx1*ly1*lz1*lelv
        call rone(dummy,ntot)
        call setmu_svv(dummy,vx,vy,vz)

        call axhelm_svv_fluid(auvx,vx,1,1)
        call axhelm_svv_fluid(auvy,vy,1,1)
        if(if3d)then
          call axhelm_svv_fluid(auvz,vz,1,1)
        endif
      endif

      svvx = auvx(ix,iy,iz,e) / vtrans(ix,iy,iz,e,1)
      svvy = auvy(ix,iy,iz,e) / vtrans(ix,iy,iz,e,1)
      svvz = auvz(ix,iy,iz,e) / vtrans(ix,iy,iz,e,1)

      return
      end

