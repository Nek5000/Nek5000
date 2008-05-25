c-----------------------------------------------------------------------
      subroutine particle_advect_std(x,xs,npart)
      include 'SIZE'
      include 'TOTAL'
c
      real x(ndim,npart),xs(ndim,2,npart)
      common /scruz/ vxl(lx1*ly1*lz1*lelt,3)
      common /scrmg/ bw (lx1*ly1*lz1*lelt,2),xo(3)
      common /padvc/ xmx(3,0:2)
c
      if (istep.eq.0) then
         call rzero(xs,2*ndim*npart) 
         c0 = 1
         c1 = 0
         c2 = 0
         n=nx1*ny1*nz1*nelv
         xmx (1,0) = glmin(xm1,n)
         xmx (1,1) = glmax(xm1,n)
         xmx (1,2) = xmx (1,1) - xmx (1,0)
         xmx (2,0) = glmin(ym1,n)
         xmx (2,1) = glmax(ym1,n)
         xmx (2,2) = xmx (2,1) - xmx (2,0)
         xmx (3,0) = glmin(zm1,n)
         xmx (3,1) = glmax(zm1,n)
         xmx (3,2) = xmx (3,1) - xmx (3,0)
      elseif (istep.eq.1) then
         c0 = 1.5
         c1 = -.5
         c2 = 0
      else
         c0 = 23.
         c1 = -16.
         c2 = 5.
         c0 = c0/12.
         c1 = c1/12.
         c2 = c2/12.
      endif
c
      call mapleg(vxl(1,1),vx,bw)    ! Map to Legendre space
      call mapleg(vxl(1,2),vy,bw)
      if (if3d) call mapleg(vxl(1,3),vz,bw)
c
      do i=1,npart
         xo(1) = x(1,i)
         xo(2) = x(2,i)
         xo(3) = x(3,i)
         do k=1,ndim
            x (k,i)   = x(k,i) + dt*(c1*xs(k,1,i) + c2*xs(k,2,i))
            xs(k,2,i) = xs(k,1,i)
            call interpl_gtp(xs(k,1,i),xo,vxl(1,k))  ! xs_0
            x (k,i)   = x(k,i) + dt*c0*xs(k,1,i)
            if (x(k,i).lt.xmx(k,0)) x(k,i) = x(k,i) + xmx(k,2) !  periodicity
            if (x(k,i).gt.xmx(k,1)) x(k,i) = x(k,i) - xmx(k,2) !  only
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine interpl_gtp(ust,xst,u)
c
c     Given u() u* := u(X*)
c
c     This form assumes:
c
c        . u is stored in a Legendre basis (NOT the std. Lagrange basis)
c
c        . the geometry is globally a rectalinear tensor-product 
c          array of elements [ nelx x nely x nelz ], and that
c          gfdm_set_geom has already been called.
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'ZPER'
c
      real u(nx1*ny1*nz1,nelv),xst(3)
      integer e,eg,ex,ey,ez,e_old,nid_old
      save nid_old,e_old
      data nid_old,e_old / 0,1 /
c
      real hrc(100),hsc(100),wk(100)
c
      real xold(3) , rold(3)
      save xold    , rold
      data xold    , rold / 6*1.e19 /
c
      ust = 0.  ! default for not on this processor
c
      if ( xst(1).eq.xold(1)   .and. 
     $     xst(2).eq.xold(2)   .and.
     $     xst(3).eq.xold(3) ) then
         if (nid.eq.nid_old) call evalscl(ust,u(1,e_old),rold)
         return
      endif
c
c     If we get here, then we have a new point.
c
      xold(1) = xst(1)
      xold(2) = xst(2)
      xold(3) = xst(3)
c
      ex = interval_s(xst(1),xgtp,nelx)
      ey = interval_s(xst(2),ygtp,nely)
      ez = interval_s(xst(3),zgtp,nelz)
      if (.not.if3d) ez = 1
      eg = ex + nelx*(ey-1) + nely*nelx*(ez-1)
c
      e_old   = gllel (eg)
      nid_old = gllnid(eg)
c
      if (nid_old.eq.nid) then
         dx = xm1(nx1,1,1,e_old) - xm1(1,1,1,e_old)
         dy = ym1(1,ny1,1,e_old) - ym1(1,1,1,e_old)
         dz = zm1(1,1,nz1,e_old) - zm1(1,1,1,e_old)
         if (.not.if3d) dz = 1
c
         rold(1)  = -1 + 2*(xst(1)-xm1(1,1,1,e_old))/dx
         rold(2)  = -1 + 2*(xst(2)-ym1(1,1,1,e_old))/dy
         rold(3)  = -1 + 2*(xst(3)-zm1(1,1,1,e_old))/dz
c
         call evalscl(ust,u(1,e_old),rold) ! evaluate w/ legendre func.
c
      endif
c
      return
      end
c-----------------------------------------------------------------------
      function interval_s (x,xx,n)
c
c     Return the interval containing x,  assuming xx_i increases monotonically
c
      real xx(0:n)
c
      i = 0
c
      do i=0,n
         if (x.le.xx(i)) then
           interval_s = i
           return
         endif
      enddo
      interval_s = n+1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine evalscl( u0 , u , r0 )
C
C     Evaluate a scalar, SCAL, at position RRL and return the result in X0.
C
      include 'SIZE'
      include 'INPUT'
c
      real u(lx1,ly1,lz1),r0(3)
c
      common  /ceval/ hr(lx1),hs(ly1),ht(lz1)
     $              , hh(ly1),hhh(ly1,lz1)
c
      real rold(3)
      save rold
      data rold / 3*1.e9 /
c
      u0   = 0
      ijk  = 0
      nxyz = nx1*ny1*nz1
c
      if ( r0(1).ne.rold(1)   .or. 
     $     r0(2).ne.rold(2)   .or.
     $     r0(3).ne.rold(3) ) then
         rold(1) = r0(1)
         rold(2) = r0(2)
         rold(3) = r0(3)
         call legendre_poly(hr,r0(1),nx1-1)
         call legendre_poly(hs,r0(2),ny1-1)
         if (if3d) call legendre_poly(ht,r0(3),nz1-1)
      endif
c
      if (if3d) then
         do k=1,nz1
         do j=1,ny1
            hhh(j,k) = 0
            do i=1,nx1
               hhh(j,k) = hhh(j,k) + hr(i)*u(i,j,k)
            enddo
         enddo
         enddo
c
         do k=1,nz1
            hh(k) = 0
            do j=1,ny1
               hh(k) = hh(k) + hs(j)*hhh(j,k)
            enddo
         enddo
c
         do k=1,nz1
            u0 = u0 + ht(k)*hh(k)
         enddo
      else
         do j=1,ny1
            hh(j) = 0
            do i=1,nx1
               hh(j) = hh(j) + hr(i)*u(i,j,1)
            enddo
         enddo
c
         do j=1,ny1
            u0 = u0 + hs(j)*hh(j)
         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mapleg(uh,u,w)
c
c     Convert u to Legendre space
c
      include 'SIZE'
      include 'TOTAL'
c
      real uh(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt),w(1)
c
      common /errcmn/ Lj(lx1*lx1),Ljt(lx1*lx1)
      real Lj,Ljt
c
      integer e
c
      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) then
         icalld = 1
         call build_legend_transform(Lj,Ljt,zgm1,nx1)
      endif
c
c
      do e=1,nelt   !  Go to Legendre space
         call tensr3(uh(1,e),nx1,u(1,e),nx1,Lj,Ljt,Ljt,w) 
      enddo
c
      return
      end
c-----------------------------------------------------------------------
