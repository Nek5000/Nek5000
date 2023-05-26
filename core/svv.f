c---------------------------------------------------------------------
      subroutine setdefault_svv
      include 'SIZE'
      include 'SVV'

      do i=1,ldimt
        svv_c0(i) = 0.1
        svvcut(i) = (lx1-1.0)/2.0
        ifnlsvv(i) = .false.

        svv_k0(i) = 1.5
        svv_k1(i) = 1.5
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine modifyDer(flag)
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
      subroutine setmu_svv(phi)
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
         call cmult(svvf,svv_c0(ifield-1)/(lx1-1.0),ntot)
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
      if(ifnlsvv(ifield-1))then
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

      call modifyDer(-1)
      call gradm1(gradx,grady,gradz,phi)
      call modifyDer(1)
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
           cnl = (lx1-1.0)**svv_k1(ifield-1)
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
         threshold = 10.0**((lx1-1.)/svv_k0(ifield-1))      
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
      include 'SIZE'
      include 'TOTAL'
      include 'SVV'
c
      real q(lx1,lx1)
      real a,k
      real n
c
      n = lx1-1.0

      call ident(q,lx1)
      
      do i=1,lx1
         k = i-1
         q(i,i) = (k/n)**(svvcut(ifield-1)/2.0)
      enddo
      
      if(nid .eq. 0)then
         write(*,*)"---------Modal parameters for SVV filter---------"
         write(*,*)(q(i,i),i=1,lx1)
         write(*,*)"-------------------------------------------------"
      endif

      
      return
      end
c---------------------------------------------------------------------
