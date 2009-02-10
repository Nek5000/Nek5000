      subroutine load_fld(string,in)
      
      include 'SIZE'
      include 'INPUT'

      character*1 string(80),fout(80),BLNK
      character*6 ext
      DATA BLNK/' '/

      call blank  (initc(1),80)

      L1=0
      DO 100 I=1,80
         IF (STRING(I).EQ.BLNK) GOTO 200
         L1=I
  100 CONTINUE
  200 CONTINUE
      LEN=L1

      if(param(67).eq.6) then
        write(ext,1) in
   1    format('.f', i4.4)
      else
        write(ext,2) in
   2    format('.fld', i2.2)
      endif

      call chcopy (fout,string,len)
      call chcopy (fout(len+1),ext,6)

      call chcopy (initc(1),fout,len+6)
      call setics

      return
      end

c-----------------------------------------------------------------------
      subroutine lambda2(l2)
c
c     Generate Lambda-2 vortex of Jeong & Hussein, JFM '95
c
      include 'SIZE'
      include 'TOTAL'

      real l2(lx1,ly1,lz1,1)

      parameter (lxyz=lx1*ly1*lz1)

      real gije(lxyz,3,3)
      real vv(ldim,ldim),ss(ldim,ldim),oo(ldim,ldim),lam(ldim)


      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv

      do ie=1,nelv
         ! Compute velocity gradient tensor
         call comp_gije(gije,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)

         do l=1,nxyz
            ! decompose into symm. and antisymm. part
            do j=1,ndim
            do i=1,ndim
               ss(i,j) = 0.5*(gije(l,i,j)+gije(l,j,i))
               oo(i,j) = 0.5*(gije(l,i,j)-gije(l,j,i))
            enddo
            enddo
         
            call rzero(vv,ldim*ldim)
            do j=1,ndim
            do i=1,ndim
            do k=1,ndim
               vv(i,j) = vv(i,j) + ss(i,k)*ss(k,j) + oo(i,k)*oo(k,j)
            enddo
            enddo
            enddo

c           Solve eigenvalue problemand sort 
c           eigenvalues in ascending order.
            call find_lam3(xl2,lam,vv,w,ndim,ierr)

            l2(l,1,1,ie) = lam(2)
         enddo
      enddo

      ! smooth field
      wght = 0.3 
      ncut = 1
      call filter_s(l2,wght,ncut,'vortx') 

      return
      end


c-----------------------------------------------------------------------
      subroutine find_lam3(lam2,lam,aa,w,ndim,ierr)
      real aa(3,3),lam(3),w(3,3),lam2
c
c     Use cubic eqn. to compute roots
c
      common /ecmnr/ a,b,c,d,e,f,f2,ef,df,r0,r1
      common /ecmni/ nr
      common /ecmnl/ iffout,ifdefl
      logical        iffout,ifdefl
c
c
      iffout = .false.
      ierr = 0
c
c     2D case....
c
c
      if (ndim.eq.2) then
         a = aa(1,1)
         b = aa(1,2)
         c = aa(2,1)
         d = aa(2,2)
         aq = 1.
         bq = -(a+d)
         cq = a*d-c*b
c
         call quadratic(x1,x2,aq,bq,cq,ierr)
c 
         lam2 = min(x1,x2)
c
         return
      endif
c
c
c
c
c     Else ...  3D case....
c
c                                    a d e
c     Get symmetric 3x3 matrix       d b f
c                                    e f c
c
      a = aa(1,1)
      b = aa(2,2)
      c = aa(3,3)
      d = 0.5*(aa(1,2)+aa(2,1))
      e = 0.5*(aa(1,3)+aa(3,1))
      f = 0.5*(aa(2,3)+aa(3,2))
      ef = e*f
      df = d*f
      f2 = f*f
c
c
c     Use cubic eqn. to compute roots
c
c     ax = a-x
c     bx = b-x
c     cx = c-x
c     y = ax*(bx*cx-f2) - d*(d*cx-ef) + e*(df-e*bx)
c
      a1 = -(a+b+c)
      a2 =  (a*b+b*c+a*c) - (d*d+e*e+f*f)
      a3 =  a*f*f + b*e*e + c*d*d - a*b*c - 2*d*e*f
c
      call cubic  (lam,a1,a2,a3,ierr)
      call sort (lam,w,3)
      lam2 = lam(2)
c
      return
      end


c-----------------------------------------------------------------------
      subroutine quadratic(x1,x2,a,b,c,ierr)
c
c     Stable routine for computation of real roots of quadratic
c
      ierr = 0
      x1 = 0.
      x2 = 0.
c
      if (a.eq.0.) then
         if (b.eq.0) then
            if (c.ne.0) then
c              write(6,10) x1,x2,a,b,c
               ierr = 1
            endif
            return
         endif
         ierr = 2
         x1 = -c/b
c        write(6,11) x1,a,b,c
         return
      endif
c
      d = b*b - 4.*a*c
      if (d.lt.0) then
         ierr = 1
c        write(6,12) a,b,c,d
         return
      endif
      if (d.gt.0) d = sqrt(d)
c
      if (b.gt.0) then
         x1 = -2.*c / ( d+b )
         x2 = -( d+b ) / (2.*a)
      else
         x1 =  ( d-b ) / (2.*a)
         x2 = -2.*c / ( d-b )
      endif
c
   10 format('ERROR: Both a & b zero in routine quadratic NO ROOTS.'
     $      ,1p5e12.4)
   11 format('ERROR: a = 0 in routine quadratic, only one root.'
     $      ,1p5e12.4)
   12 format('ERROR: negative discriminate in routine quadratic.'
     $      ,1p5e12.4)
c
      return
      end


c-----------------------------------------------------------------------
      subroutine cubic(xo,ai1,ai2,ai3,ierr)
      real xo(3),ai1,ai2,ai3
      complex*16 x(3),a1,a2,a3,q,r,d,arg,t1,t2,t3,theta,sq,a13
c
c     Compute real solutions to cubic root eqn. (Num. Rec. v. 1, p. 146)
c     pff/Sang-Wook Lee  Jan 19 , 2004
c
c     Assumption is that all x's are *real*
c
      real*8 twopi
      save   twopi
      data   twopi /6.283185307179586476925286766/
c
      ierr = 0
c
      zero = 0.
      a1   = cmplx(ai1,zero)
      a2   = cmplx(ai2,zero)
      a3   = cmplx(ai3,zero)
c
      q = (a1*a1 - 3*a2)/9.
      if (q.eq.0) goto 999
c
      r = (2*a1*a1*a1 - 9*a1*a2 + 27*a3)/54.
c
      d = q*q*q - r*r
c
c     if (d.lt.0) goto 999
c
      arg   = q*q*q
      arg   = sqrt(arg)
      arg   = r/arg
c
      if (abs(arg).gt.1) goto 999
      theta = acos(abs(arg))
c
      t1    = theta / 3.
      t2    = (theta + twopi) / 3.
      t3    = (theta + 2.*twopi) / 3.
c
      sq  = -2.*sqrt(q)
      a13 = a1/3.
      x(1) = sq*cos(t1) - a13
      x(2) = sq*cos(t2) - a13
      x(3) = sq*cos(t3) - a13
c
      xo(1) = real(x(1))
      xo(2) = real(x(2))
      xo(3) = real(x(3))
c
      return
c
  999 continue   ! failed
      ierr = 1
      call rzero(x,3)

      return
      end


c-----------------------------------------------------------------------
      subroutine comp_gije(gije,u,v,w,e)
c
c                                         du_i
c     Compute the gradient tensor G_ij := ----  ,  for element e
c                                         du_j
c
      include 'SIZE'
      include 'TOTAL'

      real gije(lx1*ly1*lz1,ndim,ndim)
      real u   (lx1*ly1*lz1)
      real v   (lx1*ly1*lz1)
      real w   (lx1*ly1*lz1)

      real ur  (lx1*ly1*lz1)
      real us  (lx1*ly1*lz1)
      real ut  (lx1*ly1*lz1)

      integer e

      common /dudxyj/ jacmi(lx1*ly1*lz1,lelt)
      real jacmi,j

      integer icalld
      save    icalld
      data    icalld /-9/

      if (istep.ne.icalld) then
         n = nx1*ny1*nz1*nelt
         call invers2(jacmi,jacm1,n)
         icalld=istep
      endif

      n    = nx1-1      ! Polynomial degree
      nxyz = nx1*ny1*nz1

      if (if3d) then     ! 3D CASE

        do k=1,3
          if (k.eq.1) call local_grad3(ur,us,ut,u,n,1,dxm1,dxtm1)
          if (k.eq.2) call local_grad3(ur,us,ut,v,n,1,dxm1,dxtm1)
          if (k.eq.3) call local_grad3(ur,us,ut,w,n,1,dxm1,dxtm1)

          do i=1,nxyz
            dj = jacmi(i,e)

            ! d/dx
            gije(i,k,1) = dj*( 
     $      ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)+ut(i)*txm1(i,1,1,e))
            ! d/dy
            gije(i,k,2) = dj*( 
     $      ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)+ut(i)*tym1(i,1,1,e))
            ! d/dz
            gije(i,k,3) = dj*(
     $      ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e))

          enddo
        enddo

      elseif (ifaxis) then   ! AXISYMMETRIC CASE
            if(nid.eq.0) write(6,*) 
     &        'ABORT: comp_gije no axialsymmetric support for now'
            call exitt
      else              ! 2D CASE

        do k=1,2
          if (k.eq.1) call local_grad2(ur,us,u,n,1,dxm1,dxtm1)
          if (k.eq.2) call local_grad2(ur,us,v,n,1,dxm1,dxtm1)
          do i=1,nxyz
             dj = jacmi(i,e)
             ! d/dx
             gije(i,k,1)=dj*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))
             ! d/dy 
             gije(i,k,2)=dj*(ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e))
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine filter_s(scalar,wght,ncut,name5) ! filter scalar field 

      include 'SIZE'
      include 'TOTAL'

      real scalar(1)
      character*5 name5

      parameter (l1=lx1*lx1)
      real intdv(l1),intuv(l1),intdp(l1),intup(l1),intv(l1),intp(l1)
      save intdv    ,intuv    ,intdp    ,intup    ,intv    ,intp

      common /ctmp0/ intt
      common /screv/ wk1,wk2
      common /scrvh/ zgmv,wgtv,zgmp,wgtp,tmax(100)

      real intt (lx1,lx1)
      real wk1  (lx1,lx1,lx1,lelt)
      real wk2  (lx1,lx1,lx1)
      real zgmv (lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)


      integer icall
      save    icall
      data    icall /0/

      logical ifdmpflt

      imax = nid
      imax = iglmax(imax,1)
      jmax = iglmax(imax,1)
      if (icall.eq.0) call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)
      icall = 1

      call filterq(scalar,intv,nx1,nz1,wk1,wk2,intt,if3d,fmax)
      fmax = glmax(fmax,1)

      if (nid.eq.0) write(6,1) istep,fmax,name5
    1 format(i8,' sfilt:',1pe12.4,a10)

      return
      end

c-----------------------------------------------------------------------
