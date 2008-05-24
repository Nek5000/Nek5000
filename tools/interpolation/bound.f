      integer function bnd_hash_cnt(ie,n)

c----------------------------------------------------------
c     If this proc's bounding box --- given by bnd_all ---
c     is divided into n parts on each side, this function
c     returns the number of sub-boxes that a given element
c     overlaps (using its bounding box from bnd_el)
c----------------------------------------------------------

      include 'SIZE'
      include 'BND'

      integer ie,n
      integer d,ia,ib

      bnd_hash_cnt = 1
      do d=1,ndim
         call bnd_hash_range(ie,n,d,ia,ib)
         bnd_hash_cnt = bnd_hash_cnt * (ib-ia+1)
      enddo

      end
     
C=======================================================================

      subroutine bnd_hash_range(ie,n,d,ia,ib)

c----------------------------------------------------------
c     If this proc's bounding box --- given by bnd_all ---
c     is divided into n parts in dimension d, this function
c     returns the index of the first and last part,
c     in ia and ib, that a given element ie overlaps
c----------------------------------------------------------

      include 'SIZE'
      include 'BND'

      integer ie,n,d,ia,ib

      ia = int(1+bnd_fac(d)*(bnd_el(1,d,ie)-bnd_all(1,d)))
      ib = n-int(bnd_fac(d)*(bnd_all(2,d)-bnd_el(2,d,ie)))
      ia = min(max(ia,1),n)
      ib = min(max(ib,1),n)
      end

C=======================================================================

      subroutine bnd_hash_inc(a,n,ia,ib,ja,jb,ka,kb)

      integer n,ia,ib,ja,jb,ka,kb
      integer a(n,n,n)
      
      integer i,j,k

      do k=ka,kb
        do j=ja,jb
          do i=ia,ib
             a(i,j,k)=a(i,j,k)+1
          enddo
        enddo
      enddo

      end

C=======================================================================

      subroutine bnd_hash_setup

c----------------------------------------------------------
c     First calculates (via a call to bnd_calc) bounding
c     boxes for all elements and for the proc as a whole.
c
c     Then sets up a compressed table with use as follows.
c
c     The proc's bbox is subdivided into n = bnd_hn pieces
c     on each side, for a total of n**ndim sub-boxes. Each
c     sub-box (i,j,k) is assigned the single index
c        ii = ((k-1)*n+(j-1))*n+i .
c     If a point is in sub-box ii, then the point could
c     possibly be in one of the elements
c        bnd_hash(bnd_hash(ii))
c        bnd_hash(bnd_hash(ii)+1)
c        ...
c        bnd_hash(bnd_hash(ii+1)-1)
c     according to the calculated element bboxes.
c
c     n = bnd_hn is chosen to be as large as possible
c     given the storage limit --- bnd_hash has size
c     lbndhmax (in BNDEL)
c----------------------------------------------------------

      include 'SIZE'
      include 'INPUT'
      include 'BND'

      integer bnd_hash_cnt          ! a function
      integer M
      parameter (M=lbndhmax)
      real one, invdim
      parameter (one=1.0,invdim=one/ldim)
      
      integer base_off
      integer ia,ib,ja,jb,ka,kb,i,j,k,ii,jj,kk
      integer nl, nm, nu, stg, ie

      call bnd_calc

      ! determine optimal n = bnd_hn using bisection
      nl = 1                          ! this will definitely fit
      nu = M-int(M-(M-nelv)**invdim)  ! this definitely won't
c      print *, 'nl/nu/M',nl,nu,M
      if(nu.lt.2) nu=2
      do while (nu-nl.gt.1)
         nm = (nl+nu)/2
         do i=1,ndim
            bnd_fac(i) = nm/(bnd_all(2,i)-bnd_all(1,i))
         enddo
         ! compute storage needed for n = nm
         stg = nm**ndim+1
         do ie=1,nelv
           stg=stg+bnd_hash_cnt(ie,nm)
         enddo
c         print *, 'nm/stg',nm,stg
         if(stg.le.M) nl = nm ! fits
         if(stg.gt.M) nu = nm ! doesn't fit
      enddo
      bnd_hn = nl
      base_off = nl**ndim+1
      do i=1,ndim
         bnd_fac(i) = nl/(bnd_all(2,i)-bnd_all(1,i))
      enddo

      ! now determine the length of each sub-box's list
      call izero(bnd_hash,M)
      do ie=1,nelv
         call bnd_hash_range(ie,nl,1,ia,ib)
         call bnd_hash_range(ie,nl,2,ja,jb)
         ka = 1
         kb = 1
         if(if3d) call bnd_hash_range(ie,nl,3,ka,kb)
         call bnd_hash_inc(bnd_hash(2),nl,ia,ib,ja,jb,ka,kb)
      enddo
      !now bnd_hash(ii+1) = len of list for sub-box ii
      !the +1 in bnd_hash(ii+1) is acheived by passing
      !  bnd_hash(2) instead of bnd_hash(1) to the
      !  increment routines
      !change to offsets
      bnd_hash(1) = base_off+1
      do i=2,base_off
         bnd_hash(i) = bnd_hash(i-1) + bnd_hash(i)
      enddo
c      print *, 'lastoff', bnd_hash(base_off)
c      call exitt
      !fill the lists
      do ie=nelv,1,-1
c        write(6,*) ie,'bnd_hash',ia,ib,ka,kb
         call bnd_hash_range(ie,nl,1,ia,ib)
         call bnd_hash_range(ie,nl,2,ja,jb)
         ka = 1
         kb = 1
         if(if3d) call bnd_hash_range(ie,nl,3,ka,kb)
         do k=ka,kb
         do j=ja,jb
         do i=ia,ib
            ii=((k-1)*bnd_hn+(j-1))*bnd_hn+i
            jj=bnd_hash(ii)
            len=bnd_hash(ii+1)-jj
            kk=bnd_hash(jj)
            bnd_hash(jj+len-kk-1)=ie
            if(kk.lt.len-1) bnd_hash(jj)=kk+1
         enddo
         enddo
         enddo
      enddo
      end

C=======================================================================

      subroutine bnd_calc

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'BND'

      integer ie,d

      call bnd_setup
      if(nelv.eq.0) return
      if(.not.if3d) call bnd_el2d(xm1,ym1,bnd_el(1,1,1))
      if(     if3d) call bnd_el3d(xm1,ym1,zm1,bnd_el(1,1,1))
      call copy(bnd_all,bnd_el,2*ndim)
      nxyz = nx1*ny1*nz1
      do ie=2,nelv
         !if (mod(ie,100).eq.0) write(6,*) ie,' bnd_calc',ioff
         ioff = nxyz*(ie-1) + 1
         if(.not.if3d) call bnd_el2d(xm1(ioff)
     $      ,ym1(ioff),bnd_el(1,1,ie))
         if(     if3d) call bnd_el3d(xm1(ioff)
     $      ,ym1(ioff),zm1(ioff),bnd_el(1,1,ie))
         do d=1,ndim
c            print *, ie,bnd_el(1,d,ie),bnd_el(2,d,ie)
            bnd_all(1,d)=min(bnd_all(1,d),bnd_el(1,d,ie))
            bnd_all(2,d)=max(bnd_all(2,d),bnd_el(2,d,ie))
         enddo
c        xmx = max(xmx,glmax(zm1(ioff),nxyz))
c        write(6,*) ie,ioff,xmx,' xmx'
      enddo

      do d=1,ndim
         print *, 'tbound',d,bnd_all(1,d),bnd_all(2,d)
      enddo
c     call exitt

      end

C=======================================================================

      subroutine bnd_setup

      include 'SIZE'
      include 'BND'

      real            Jm               , Dm
      common /bndwk1/ w(lx1)
      common /bndwka/ Jm(2*lbndn-1,lx1), Dm(2*lbndn-1,lx1)

      real pi
      parameter (pi=3.14159265358979323846264338)

      integer i,j,n,m
      
      n=nx1
      m=lbndn
c     bnd_m = m

      bnd_m = min(m,2*n)
      m=bnd_m

c     bnd_z are the nx1   GLL nodes and
      call zwgll(bnd_z,w,n)
      call fdwquick_init(bnd_z,n,bnd_fdw)

c     bnd_h are the bnd_m Chebyshev nodes

      bnd_h(1) = -1.0
      bnd_h(m) =  1.0
      do i=2,m-1
         bnd_h(i) = cos((m-i)*pi/(m-1))
      enddo
      
c     Helper routine finds uv, ov
c        uv(:,j) and ov(:,j) are nodal values (bnd_h nodes)
c          for lower and upper piecewise linear bounds on
c          the jth GLL polynomial Lagrangian basis function
c     Also finds uvp, uvn, ovp, ovn
c        which are the positive and negative parts of uv, ov

      call bnd_setup_hlpr(bnd_z,n,bnd_h,m,Jm,Dm)
      call gll2leg(n,bnd_z,w,2,Qm)

      end
            
C=======================================================================

      subroutine bnd_setup_hlpr(z,n,h,m,Jm,Dm)

      integer n,m
      real z(n), h(m)
      
      include 'SIZE'
      include 'BND'

      real Jm(2*m-1,lx1), Dm(2*m-1,lx1)
      common /bndwk2/ q(2*lbndn-1), work(lx1,5)
      
      integer i,ii,j,mm
      real c0,c1
      
      mm = 2*m-1
      do i=1,m-1
         q(2*i-1) = h(i)
         q(2*i) = 0.5*(h(i)+h(i+1))
      enddo
      q(mm) = h(m)
      
      call fdwquick01(z,n,q,mm,Jm,Dm,bnd_fdw,work)

c     find piecewise linear uv(:,j), ov(:,j) such that
c     uv(:,j) < pi(j) < ov(:,j)
c     where pi(j) is the jth Lagrangian basis function

      do j=1,n
         uv(1,j) = Jm(1,j)
         ov(1,j) = Jm(1,j)
         uv(m,j) = Jm(mm,j)
         ov(m,j) = Jm(mm,j)
         do i=2,m-1
            ii = 2*i-1
            c0 = Jm(ii-1,j) + (q(ii)-q(ii-1))*Dm(ii-1,j)
            c1 = Jm(ii+1,j) + (q(ii)-q(ii+1))*Dm(ii+1,j)
            uv(i,j) = min(Jm(ii,j),c0,c1)
            ov(i,j) = max(Jm(ii,j),c0,c1)
         enddo
      enddo
      
c     split uv, ov into positive and negative parts      

      do j=1,n
         do i=1,m
            uvp(i,j)=0.0
            ovp(i,j)=0.0
            uvn(i,j)=0.0
            ovn(i,j)=0.0
            if(uv(i,j).gt.0) uvp(i,j)=uv(i,j)
            if(uv(i,j).lt.0) uvn(i,j)=uv(i,j)
            if(ov(i,j).gt.0) ovp(i,j)=ov(i,j)
            if(ov(i,j).lt.0) ovn(i,j)=ov(i,j)
         enddo
      enddo      

      end

C=======================================================================

      subroutine bnd_lines1d(u,su,a,sa,b,sb)

      include 'SIZE'
      include 'BND'

      integer su,sa,sb   ! strides
      real u(su,nx1), a(sa,bnd_m), b(sb,bnd_m)
      
      real a0, a1, w
      integer i
      
c     Compute first two coefficients in Legendre expansion
      a0 = 0
      a1 = 0
      do i=1,nx1
         a0 = a0 + Qm(1,i)*u(1,i)
         a1 = a1 + Qm(2,i)*u(1,i)
      enddo
c     a0 and a1 can be set to anything, but the bounds will
c        be better when u(x) - a0-a1*x is smaller
c     Legendre coefficients seem like a good choice
      do i=1,bnd_m
         a(1,i) = a0 + a1*bnd_h(i)
         b(1,i) = a(1,i)
      enddo
      
      do i=1,nx1
         w = u(1,i) - a0 - a1*bnd_z(i)
         if(w.ge.0) then
            do j=1,bnd_m
               a(1,j) = a(1,j) + w*uv(j,i)
               b(1,j) = b(1,j) + w*ov(j,i)
            enddo
         else
            do j=1,bnd_m
               a(1,j) = a(1,j) + w*ov(j,i)
               b(1,j) = b(1,j) + w*uv(j,i)
            enddo
         endif
      enddo
      
      end

C=======================================================================

      subroutine bnd_lines1dvar(uu,suu,ou,sou,a,sa,b,sb)
      
      include 'SIZE'
      include 'BND'

      integer suu,sou,sa,sb             ! strides
      real uu(suu,nx1), ou(sou,nx1)     ! bounds on Lag. comps
      real a(sa,bnd_m), b(sb,bnd_m)     ! piecewise lin. bounds
      
      real a0, a1, temp, uw, ow
      integer i

      real half
      parameter (half=0.5)
      
c     Compute first two coefficients in Legendre expansion
c        of "average" polynomial
      a0 = 0.0
      a1 = 0.0
      do i=1,nx1
         temp = half*(uu(1,i)+ou(1,i))
         a0 = a0 + Qm(1,i)*temp
         a1 = a1 + Qm(2,i)*temp
      enddo
c     a0 and a1 can be set to anything, but the bounds will
c        be better when {uu(x),ou(x)} - a0-a1*x is smaller
c     Legendre coefficients seem like a good choice
      do i=1,bnd_m
         a(1,i) = a0 + a1*bnd_h(i)
         b(1,i) = a(1,i)
      enddo
      
      do i=1,nx1
         temp = a0+a1*bnd_z(i)
         uw = uu(1,i) - temp
         ow = ou(1,i) - temp
         if(uw.ge.0.0) then     !  0 <= uw <= ow
            do j=1,bnd_m
               a(1,j) = a(1,j) + uw*uvp(j,i) + ow*uvn(j,i)
               b(1,j) = b(1,j) + ow*ovp(j,i) + uw*ovn(j,i)
            enddo
         elseif(ow.le.0.0) then ! uw <= ow <= 0
            do j=1,bnd_m
               a(1,j) = a(1,j) + uw*ovp(j,i) + ow*ovn(j,i)
               b(1,j) = b(1,j) + ow*uvp(j,i) + uw*uvn(j,i)
            enddo
         else                 ! uw <  0  <  ow
            do j=1,bnd_m
               a(1,j) = a(1,j) + uw*ovp(j,i) + ow*uvn(j,i)
               b(1,j) = b(1,j) + ow*ovp(j,i) + uw*uvn(j,i)
            enddo
         endif
      enddo
      
      end

C=======================================================================

      subroutine bnd_1d(u,su,a,b)

      include 'SIZE'
      include 'BND'

      integer su                         ! stride
      real u(su,nx1), a, b
      
      common /bndwk2/ av(lbndn), bv(lbndn)
      real av,bv
      
      integer i
      
      call bnd_lines1d(u,su,av,1,bv,1)
      a = av(1)
      b = bv(1)
      do i=2,bnd_m
         a = min(a,av(i))
         b = max(b,bv(i))
      enddo
      end

C=======================================================================

      subroutine bnd_2d(u,su1,su2,a,b)

      include 'SIZE'
      include 'BND'

      integer su1,su2   ! strides
      real u(su1,nx1,su2,nx1), a, b
      
      real av,bv,uu,ou
      common /bndwk2/ av(lx1,lbndn), bv(lx1,lbndn),
     $                uu(lbndn,lbndn), ou(lbndn,lbndn)
      
      integer i,j
      do i=1,nx1
         call bnd_lines1d(u(1,1,1,i),su1,av(i,1),lx1,bv(i,1),lx1)
      enddo
      do i=1,bnd_m
         call bnd_lines1dvar(av(1,i),1,bv(1,i),1,uu(1,i),1,ou(1,i),1)
      enddo
      
      a = uu(1,1)
      b = ou(1,1)
      do j=1,bnd_m
         do i=1,bnd_m
            a = min(a,uu(i,j))
            b = max(b,ou(i,j))
         enddo
      enddo
      end

C=======================================================================

      subroutine bnd_coord2d(x,ca,cb)

      include 'SIZE'

      real x(nx1,nx1), ca, cb

      real a, b
      common /bndwk1/ a(4), b(4)
      real dist

      call bnd_1d(x(  1,  1),  1,a(1),b(1))
      call bnd_1d(x(  1,nx1),  1,a(2),b(2))
      call bnd_1d(x(  1,  1),nx1,a(3),b(3))
      call bnd_1d(x(nx1,  1),nx1,a(4),b(4))
      ca = min(a(1),a(2),a(3),a(4))
      cb = max(b(1),b(2),b(3),b(4))
      dist = 1.e-6*(cb-ca)
      ca = ca - dist                 ! failsafe
      cb = cb + dist                 ! failsafe
      
      end
C=======================================================================
      
      subroutine bnd_el2d(x,y,bound)

      real x(1), y(1), bound(2,2)

      call bnd_coord2d(x,bound(1,1),bound(2,1))
      call bnd_coord2d(y,bound(1,2),bound(2,2))

      end

C=======================================================================

      subroutine bnd_coord3d(x,ca,cb)

      include 'SIZE'

      real x(nx1,nx1,nx1), ca, cb

      real a, b
      common /bndwk1/ a(6), b(6)      
      real dist
      
      call bnd_2d(x(  1,  1,  1),  1,  1,a(1),b(1))
      call bnd_2d(x(  1,  1,nx1),  1,  1,a(2),b(2))
      call bnd_2d(x(  1,  1,  1),  1,nx1,a(3),b(3))
      call bnd_2d(x(  1,nx1,  1),  1,nx1,a(4),b(4))
      call bnd_2d(x(  1,  1,  1),nx1,  1,a(5),b(5))
      call bnd_2d(x(nx1,  1,  1),nx1,  1,a(6),b(6))
      ca = min(a(1),a(2),a(3),a(4),a(5),a(6))
      cb = max(b(1),b(2),b(3),b(4),b(5),b(6))
      dist = 1.0e-6*(cb-ca)
      ca = ca - dist   ! failsafe
      cb = cb + dist   ! failsafe

      end

C=======================================================================

      subroutine bnd_el3d(x,y,z,bound)

      real x(1), y(1), z(1), bound(2,3)

      call bnd_coord3d(x,bound(1,1),bound(2,1))
      call bnd_coord3d(y,bound(1,2),bound(2,2))
      call bnd_coord3d(z,bound(1,3),bound(2,3))

      end

C=======================================================================

      subroutine gll2leg(n,z,w,m,Q)

c----------------------------------------------------------
c     GLL Lagrangian to Legendre Change-of-basis Matrix
c
c     Given: GLL nodes and weights z(1:n), w(1:n)
c     Computes: m rows of change-of-basis matrix Q(0:m-1,n)
c     Assumes: 2 <= m <= n
c             n
c            __                     th 
c     a   =  >   Q   u      is the i   Legendre component
c      i     --   ij  j       given the GLL Lagrangian
c            j=1              components u
c                                         j
c----------------------------------------------------------

      implicit none

      integer n,m
      real z(n), w(n), Q(0:m-1,n)
      
      real P(0:1)
      integer j,i,l,lp1,mm
      
      real one,half
      parameter (one=1.0, half=0.5)
      
      mm = min(m-2,n-3)
      do j=1,n
         P(0)=1
         P(1)=z(j)
         Q(0,j) = half*w(j)
         Q(1,j) = 3.0*half*w(j)*z(j)
         do i=1,mm
            l = iand(i,1)
            lp1 = ieor(l,1)
            P(lp1) = ( (2.0*i+1.0)*z(j)*P(l) - i*P(lp1) )/(i+1.0)
            Q(i+1,j) = half*(2*i+3)*w(j)*P(lp1)
         enddo
      enddo
      if(m.eq.n) then
         do j=1,n
            Q(n-1,j)=0
            if(j.eq.n) Q(n-1,j)=1
            do i=0,n-2
               Q(n-1,j) = Q(n-1,j)-Q(i,j)
            enddo
         enddo
      endif
      
      end
