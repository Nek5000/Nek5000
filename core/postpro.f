      subroutine load_fld(string)

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

      call chcopy (initc(1),string,len)
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
      real vv(ldim,ldim),ss(ldim,ldim),oo(ldim,ldim),w(ldim,ldim)
      real lam(ldim)

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
            call find_lam3(lam,vv,w,ndim,ierr)

            l2(l,1,1,ie) = lam(2)
         enddo
      enddo

      ! smooth field
      wght = 0.5 
      ncut = 1
      call filter_s0(l2,wght,ncut,'vortx') 

      return
      end
c-----------------------------------------------------------------------
      subroutine find_lam3(lam,aa,w,ndim,ierr)
      real aa(ndim,ndim),lam(ndim),w(ndim,ndim),lam2
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
         lam(1) = min(x1,x2)
         lam(2) = max(x1,x2)
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
      call sort   (lam,w,3)
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

      real gije(lx1*ly1*lz1,ldim,ldim)
      real u   (lx1*ly1*lz1)
      real v   (lx1*ly1*lz1)
      real w   (lx1*ly1*lz1)

      real ur  (lx1*ly1*lz1)
      real us  (lx1*ly1*lz1)
      real ut  (lx1*ly1*lz1)

      integer e

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
      subroutine filter_s1(scalar,tf,nx,nel) ! filter scalar field 

      include 'SIZE'

      parameter(lxyz=lx1*ly1*lz1) 
      real scalar(lxyz,1)
      real fh(nx*nx),fht(nx*nx),tf(nx)

      real w1(lxyz,lelt)

c     Build 1D-filter based on the transfer function (tf)
      call build_1d_filt(fh,fht,tf,nx,nid)

c     Filter scalar
      call copy(w1,scalar,lxyz*nel)
      do ie=1,nel
         call tens3d1(scalar(1,ie),w1(1,ie),fh,fht,nx1,nx1)  ! fh x fh x fh x scalar
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine filter_s0(scalar,wght,ncut,name5) ! filter scalar field 

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

c      if (icall.eq.0) call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)
      call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)

      icall = 1

      call filterq(scalar,intv,nx1,nz1,wk1,wk2,intt,if3d,fmax)
      fmax = glmax(fmax,1)

      if (nid.eq.0) write(6,1) istep,fmax,name5
    1 format(i8,' sfilt:',1pe12.4,a10)

      return
      end
c-----------------------------------------------------------------------
      subroutine intpts_setup(bb_t)
c IN:
c bb_t ... bounding box tolerance (relative to the element size)

      INCLUDE 'SIZE'
      INCLUDE 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
      common /intp/   icrh,ipth,noff,lxyz,idim

      lxyz = lx1*ly1*lz1
      noff = lxyz*lelt 
      idim = ndim
c
      call crystal_new(icrh,nekcomm,npp)
      call findpts_new(ipth,icrh,ndim,xm1,ym1,zm1,nx1,ny1,nz1,nelt,bb_t)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine intpts(field,nfld,iTupleList,mi,rTupleList,mr,n,nmax)
c IN:
c field   ... field(s) to interpolate
c nfld    ... number of fields
c n       ... local number of interpolation points
c nmax    ... maximum number of local points
c             (nmax=[n,gsum(n)])
c
c IN/OUT:
c Each interpolation point i (i=1,...,n) comes with a tuple list 
c TYPE: integer iTupleList(4,nmax), real rTupleList(1+2*ndim+nfld)
c
c iTupleList  ... integer tuple list
c    output   vi(1,i) = processor number (0 to np-1)   
c    output   vi(2,i) = local element number (1 to nelt)
c    output   vi(3,i) = return code (-1, 0, 1)
c    output   vi(4,i) = local point id (only internally used)
c rTupleList  ... real tuple list  
c    output   vr(1,i)      = distance (from located point to given point)
c    input    vr(2,i)      = x  
c    input    vr(3,i)      = y  
c    input    vr(4,i)      = z  (only when ndim=3)
c    output   vr(ndim+2,i) = r   
c    output   vr(ndim+3,i) = s   
c    output   vr(ndim+4,i) = t  (only when ndim=3)
c    output   vr(1+2*ndim+ifld,i) = interpolated field value
c                                   (ifld=1,nfld)
c
      real    field (*)
      integer iTupleList (mi,1)
      real    rTupleList (mr,1)

      common /intp/ icrh,ipth,noff,lxyz,idim

      if(mi.lt.4 .or. mr.lt.1+2*idim+nfld) then
        write(6,*) 'ABORT: intpts() invalid tuple size mi/mir', mi, mr
        call exitt
      endif

      ! locate points (iel,iproc,r,s,t)
      iguess = 0 ! no guess
      call findpts(ipth,n,iTupleList,mi,rTupleList,mr,iguess)

      do in=1,n
         ! store local id to preserve ordering
         iTupleList(4,in) = in 
         ! check return code 
         if(iTupleList(3,in).eq.1) then
           dist = rTupleList(1,in)
           write(6,'(A,4E15.7)') 
     &       'WARNING: point on boundary or outside the mesh xy[z]d: ',
     &       (rTupleList(1+k,in),k=1,idim),dist
         elseif(iTupleList(3,in).eq.-1) then
           write(6,'(A,3E15.7)') 
     &       'WARNING: point not within mesh xy[z]: !',
     &       (rTupleList(1+k,in),k=1,idim)
         endif
      enddo

      ! transfer point (tuple list) to the target proc that owns it
      nin = n 
      call findpts_transfer(ipth,nin,nmax,iTupleList,mi,rTupleList,mr)

      if(nin.eq.nmax+1) then ! check for error condition
        write(6,*) 'ABORT: intpts() more local points than nmax.'
        call exitt
      endif

      ! every point is local now, do interpolation ...
      do ifld = 1,nfld
         do in   = 1,nin
            iel  = iTupleList(2,in)
            ioff = (ifld-1)*noff + (iel-1)*lxyz 
            r    = rTupleList(idim+2,in)
            s    = rTupleList(idim+3,in)
            t    = rTupleList(idim+4,in) 
            call findpts_weights(ipth,r,s,t)
            call findpts_eval(ipth,rTupleList(1+2*idim+ifld,in),
     &                        field(ioff+1))
         enddo
      enddo

      ! tranfer points back to the source proc
      ! NOTE: iTupleList(1,:) has the source proc after findpts_transfer()
      call findpts_transfer(ipth,nin,nmax,iTupleList,mi,rTupleList,mr)

      if(nin.eq.nmax+1) then
        write(6,*) 'ABORT: intpts() more incoming points than nmax.'
        call exitt
      endif

      ! Restore initial tuple list ordering
      ikey = 4 ! sort index is local point id  
      call ftuple_list_sort(nin,ikey,iTupleList,mi,rTupleList,mr)


      return
      end
c-----------------------------------------------------------------------
      subroutine intpts_done()

      common /intp/ icrh,ipth,noff,lxyz,idim

      call findpts_done(ipth)
      call crystal_done(icrh)

      return
      end
c-----------------------------------------------------------------------
      subroutine tens3d1(v,u,f,ft,nv,nu)  ! v = F x F x F x u

c     Note: this routine assumes that nx1=ny1=nz1
c
      include 'SIZE'
      include 'INPUT'

      parameter (lw=4*lx1*lx1*lz1)
      common /ctensor/ w1(lw),w2(lw)

      real v(nv,nv,nv),u(nu,nu,nu)
      real f(1),ft(1)

      if (nu*nu*nv.gt.lw) then
         write(6,*) nid,nu,nv,lw,' ERROR in tens3d1. Increase lw.'
         call exitt
      endif

      if (if3d) then
         nuv = nu*nv
         nvv = nv*nv
         call mxm(f,nv,u,nu,w1,nu*nu)
         k=1
         l=1
         do iz=1,nu
            call mxm(w1(k),nv,ft,nu,w2(l),nv)
            k=k+nuv
            l=l+nvv
         enddo
         call mxm(w2,nvv,ft,nu,v,nv)
      else
         call mxm(f ,nv,u,nu,w1,nu)
         call mxm(w1,nv,ft,nu,v,nv)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine build_1d_filt(fh,fht,trnsfr,nx,nid)
c
c     This routing builds a 1D filter with transfer function diag()
c
c     Here, nx = number of points
c
      real fh(nx,nx),fht(nx,nx),trnsfr(nx)
c
      parameter (lm=40)
      parameter (lm2=lm*lm)
      common /cfiltr/ phi(lm2),pht(lm2),diag(lm2),rmult(lm),Lj(lm)
     $              , zpts(lm)
      real Lj

      common /cfilti/ indr(lm),indc(lm),ipiv(lm)
c
      if (nx.gt.lm) then
         write(6,*) 'ABORT in set_filt:',nx,lm
         call exitt
      endif

      call zwgll(zpts,rmult,nx)

      kj = 0
      n  = nx-1
      do j=1,nx
         z = zpts(j)
         call legendre_poly(Lj,z,n)
         kj = kj+1
         pht(kj) = Lj(1)
         kj = kj+1
         pht(kj) = Lj(2)
         do k=3,nx
            kj = kj+1
            pht(kj) = Lj(k)-Lj(k-2)
         enddo
      enddo
      call transpose (phi,nx,pht,nx)
      call copy      (pht,phi,nx*nx)
      call gaujordf  (pht,nx,nx,indr,indc,ipiv,ierr,rmult)

      call rzero(diag,nx*nx)
      k=1 
      do i=1,nx
         diag(k) = trnsfr(i)
         k = k+(nx+1)
      enddo

      call mxm  (diag,nx,pht,nx,fh,nx)      !          -1
      call mxm  (phi ,nx,fh,nx,pht,nx)      !     V D V

      call copy      (fh,pht,nx*nx)
      call transpose (fht,nx,fh,nx)

      do k=1,nx*nx
         pht(k) = 1.-diag(k)
      enddo
      np1 = nx+1
      if (nid.eq.0) then
         write(6,6) 'flt amp',(pht (k),k=1,nx*nx,np1)
         write(6,6) 'flt trn',(diag(k),k=1,nx*nx,np1)
   6     format(a8,16f7.4,6(/,8x,16f7.4))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mag_tensor_e(mag,aije)
c
c     Compute magnitude of tensor A_e for element e
c
c     mag(A_e) = sqrt( 0.5 (A:A) )
c
      include 'SIZE'
      REAL mag (lx1*ly1*lz1)
      REAL aije(lx1*ly1*lz1,ldim,ldim)

      nxyz = nx1*ny1*nz1

      call rzero(mag,nxyz)
 
      do 100 j=1,ndim
      do 100 i=1,ndim
      do 100 l=1,nxyz 
         mag(l) = mag(l) + 0.5*aije(l,i,j)*aije(l,i,j)
 100  continue

      call vsqrt(mag,nxyz)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_sije(gije)
c
c     Compute symmetric part of a tensor G_ij for element e
c
      include 'SIZE'
      include 'TOTAL'

      real gije(lx1*ly1*lz1,ldim,ldim)

      nxyz = nx1*ny1*nz1

      k = 1

      do j=1,ndim
      do i=k,ndim
         do l=1,nxyz
            gije(l,i,j) = 0.5*(gije(l,i,j)+gije(l,j,i))
            gije(l,j,i) = gije(l,i,j)
         enddo
      enddo
         k = k + 1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wlapl(out1,out2,out3,out4,u,h1,nel)

      include 'SIZE'
      include 'TOTAL'

      REAL           out1  (LX1,LY1,LZ1,1)
     $ ,             out2  (LX1,LY1,LZ1,1)
     $ ,             out3  (LX1,LY1,LZ1,1)
     $ ,             out4  (LX1,LY1,LZ1,1)
     $ ,             u     (LX1,LY1,LZ1,1)
     $ ,             h1    (LX1,LY1,LZ1,1)
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
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV


      NXY=NX1*NY1
      NYZ=NY1*NZ1
      NXZ=NX1*NZ1
      NXYZ=NX1*NY1*NZ1
      NTOT=NXYZ*NEL

      CALL RZERO (out1,NTOT)
      CALL RZERO (out2,NTOT)
      CALL RZERO (out3,NTOT)
      CALL RZERO (out4,NTOT)

      DO 100 IEL=1,NEL
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
C
        IF (NDIM.EQ.2) THEN
C       2-d case ...............
C       General case, speed-up for undeformed elements
C
           CALL MXM  (DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           CALL MXM  (U(1,1,1,IEL),NX1,DYTM1,NY1,DUDS,NY1)
           CALL COL3 (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3 (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)
           IF (IFDFRM(IEL)) THEN
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
           ENDIF
           CALL COL2 (TMP1,H1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,H1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           CALL MXM  (TMP2,NX1,DYM1,NY1,TM2,NY1)
           CALL ADD2 (OUT1(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (OUT2(1,1,1,IEL),TM2,NXYZ)
           CALL ADD3 (OUT4(1,1,1,IEL),TM1,TM2,NXYZ)
        ELSE
C       3-d case ...............
C       General case, speed-up for undeformed elements
C
           CALL MXM(DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           DO 10 IZ=1,NZ1
              CALL MXM(U(1,1,IZ,IEL),NX1,DYTM1,NY1,DUDS(1,1,IZ),NY1)
   10      CONTINUE
           CALL MXM     (U(1,1,1,IEL),NXY,DZTM1,NZ1,DUDT,NZ1)
           CALL COL3    (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3    (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)
           CALL COL3    (TMP3,DUDT,G3M1(1,1,1,IEL),NXYZ)
           IF (IFDFRM(IEL)) THEN
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
           ENDIF
           CALL COL2 (TMP1,H1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,H1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP3,H1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           DO 20 IZ=1,NZ1
              CALL MXM(TMP2(1,1,IZ),NX1,DYM1,NY1,TM2(1,1,IZ),NY1)
   20      CONTINUE
           CALL MXM  (TMP3,NXY,DZM1,NZ1,TM3,NZ1)
           CALL ADD2 (OUT1(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (OUT2(1,1,1,IEL),TM2,NXYZ)
           CALL ADD2 (OUT3(1,1,1,IEL),TM3,NXYZ)
           CALL ADD4 (OUT4(1,1,1,IEL),TM1,TM2,TM3,NXYZ)
        ENDIF
C
 100  CONTINUE
      

      return
      end
