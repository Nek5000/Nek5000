c---------------------------------------------------------------------
      subroutine setmu_svv(phi)
      implicit none
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

      integer nxyz,ntot
c
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv
      
      if(icalld .eq. 0)then
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
         call cmult(svvf,svv_c0/(lx1-1.0),ntot)
         icalld = 1
      endif

c     Scale with advection velocity
      call col3(svvmu,vx,vx,ntot)
      call addcol3(svvmu,vy,vy,ntot)
      if(if3d) call addcol3(svvmu,vz,vz,ntot)
      do i=1,ntot
         svvmu(i,1,1,1) = svvmu(i,1,1,1)**0.5
      enddo
      
      call col2(svvmu,svvf,ntot)
      
      return
      end
c---------------------------------------------------------------------
      subroutine modifyDer(flag)
      implicit none
c
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      integer flag

      integer icalld
      save icalld
      data icalld /0/
c            
      if(icalld .eq. 0)then
         call getconvD
         icalld = 1
      endif
      
      if(flag .eq. 1)then
         call copy(dxm1,cdxm1,lx1*lx1)
         call copy(dxtm1,cdxtm1,lx1*lx1)

         call copy(dym1,cdym1,lx1*lx1)
         call copy(dytm1,cdytm1,lx1*lx1)

         if(if3d)then
            call copy(dzm1,cdzm1,lx1*lx1)
            call copy(dztm1,cdztm1,lx1*lx1)
         endif
      elseif(flag .eq. -1)then
         call copy(dxm1,odxm1,lx1*lx1)
         call copy(dxtm1,odxtm1,lx1*lx1)

         call copy(dym1,odym1,lx1*lx1)
         call copy(dytm1,odytm1,lx1*lx1)

         if(if3d)then
            call copy(dzm1,odzm1,lx1*lx1)
            call copy(dztm1,odztm1,lx1*lx1)
         endif
      endif
      
      return
      end
c---------------------------------------------------------------------
      subroutine getconvD
      implicit none
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
      implicit none
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
      implicit none

c     Get all the 1d Legendre basis at Z for polynomial order NZ-1
      real B(NZ,NZ)
      real Z(NZ)
      real PNLEG
      integer i,j

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
      implicit none
c      
      include 'SIZE'
      include 'TOTAL'
c
      real q(lx1,lx1)
      real a,k
      real n
      real svvcut
      integer i
c
      svvcut = (lx1-1.0)/2.0
      n = lx1-1.0

      call ident(q,lx1)
      
      do i=1,lx1
         k = i-1.0
         q(i,i) = (k/n)**(svvcut/2.0)
      enddo
      
      if(nid .eq. 0)then
         write(*,*)"---------Modal parameters for SVV filter---------"
         write(*,*)(q(i,i),i=1,lx1)
         write(*,*)"-------------------------------------------------"
      endif

      
      return
      end
c---------------------------------------------------------------------
      subroutine axhelm_svv(au,u,helm1,helm2,imsh,isd)
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

      IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMSH)
      CALL RZERO (AU,NTOT)

      do 100 e=1,nel
C
        if (ifaxis) call setaxdy ( ifrzer(e) )
C
        IF (ldim.EQ.2) THEN
C
C       2-d case ...............
C
           if (iffast(e)) then
C
C          Fast 2-d mode: constant properties and undeformed element
C
           h1 = helm1(1,1,1,e)
           call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
           call mxm   (u(1,1,1,e),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,g4m1(1,1,1,e),nxyz)
           call col2  (tm2,g5m1(1,1,1,e),nxyz)
           call add3  (au(1,1,1,e),tm1,tm2,nxyz)
           call cmult (au(1,1,1,e),h1,nxyz)
C
           else
C
C
           call mxm  (cdxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           call mxm  (u(1,1,1,e),lx1,cdytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
           endif
           call col2 (tmp1,helm1(1,1,1,e),nxyz)
           call col2 (tmp2,helm1(1,1,1,e),nxyz)
           call mxm  (cdxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,cdym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)

        endif
C
        else
C
C       3-d case ...............
C
           if (iffast(e)) then
C
C          Fast 3-d mode: constant properties and undeformed element
C
           h1 = helm1(1,1,1,e)
           call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,e),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,e),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,g4m1(1,1,1,e),nxyz)
           call col2  (tm2,g5m1(1,1,1,e),nxyz)
           call col2  (tm3,g6m1(1,1,1,e),nxyz)
           call add3  (au(1,1,1,e),tm1,tm2,nxyz)
           call add2  (au(1,1,1,e),tm3,nxyz)
           call cmult (au(1,1,1,e),h1,nxyz)
C
           else
C
C
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
           call col2 (tmp1,helm1(1,1,1,e),nxyz)
           call col2 (tmp2,helm1(1,1,1,e),nxyz)
           call col2 (tmp3,helm1(1,1,1,e),nxyz)
           call mxm  (cdxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,cdym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,cdzm1,lz1,tm3,lz1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)
           call add2 (au(1,1,1,e),tm3,nxyz)
C
           endif
c
        endif
C
 100  continue
C
      if(ifield.gt.1 .and. ifsvv)call col2(au,svvmu,ntot)
         
      if (ifh2) call addcol4 (au,helm2,bm1,u,ntot)
C
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
      if (ifaxis.and.(isd.eq.2)) then
         do 200 e=1,nel
C
            if (ifrzer(e)) then
               call mxm(u  (1,1,1,e),lx1,datm1,ly1,duax,1)
               call mxm(ym1(1,1,1,e),lx1,datm1,ly1,ysm1,1)
            endif
c
            do 190 j=1,ly1
            do 190 i=1,lx1
C               if (ym1(i,j,1,e).ne.0.) then
                  if (ifrzer(e)) then
                     term1 = 0.0
                     if(j.ne.1) 
     $             term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                     term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i)
     $                       *jacm1(i,1,1,e)/ysm1(i)
                  else
                   term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                     term2 = 0.
                  endif
                  au(i,j,1,e) = au(i,j,1,e)
     $                          + helm1(i,j,1,e)*(term1+term2)
C               endif
  190       continue
  200    continue
      endif

      taxhm=taxhm+(dnekclock()-etime1)
      return
      end
