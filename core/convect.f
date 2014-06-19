c-----------------------------------------------------------------------
c
c    Stability limits:
c
c    AB3:    .7236                     w/safety (1.2):   .603
c
c    RK3:    1.73   (sqrt 3)           w/safety (1.2):   1.44
c
c    RK4:    2.828                     w/safety (1.2):   2.36
c
c    SEM Safety factor:  1.52 for N=3
c                     <  1.20 for N=16
c                     ~  1.16 for N=256
c
c-----------------------------------------------------------------------
      subroutine setup_convect(igeom)
      include 'SIZE'
      include 'TOTAL'
      logical ifnew

      common /cchar/ ct_vx(0:lorder+1) ! time for each slice in c_vx()

      if (igeom.eq.1) return
      if (param(99).lt.0) return ! no dealiasing

      if (ifchar) then

         nelc = nelv
         if (ifmhd) nelc = max(nelv,nelfld(ifldmhd))
         if (ifmhd) call exitti('no characteristics for mhd yet$',istep)

         ifnew = .true.
         if (igeom.gt.2) ifnew = .false.

         call set_conv_char(ct_vx,c_vx,vx,vy,vz,nelc,time,ifnew)

      else

         if (.not.ifpert) then
            if (ifcons) then
               call set_convect_cons (vxd,vyd,vzd,vx,vy,vz)
            else
               call set_convect_new  (vxd,vyd,vzd,vx,vy,vz)
            endif
         endif

      endif

c     write(6,*) istep,' conv',ifnew,igeom,' continu? ',time
c     read(5,*) dum

      return
      end
c-----------------------------------------------------------------------
      subroutine char_conv(p0,u,ulag,msk,c,cs,gsl)
c
c
c     Convect over last NBD steps using characteristics scheme
c
c     NOTE:  Here, we assume that ulag is stored by time-slice first,
c            then by field number (this is opposite to prior Nek5000)
c
c
      include 'SIZE'
      include 'TOTAL'
      real    p0(1),u(1),ulag(1),msk(1),c(1),cs(0:1)
      integer gsl

      common /scrns/ ct  (lxd*lyd*lzd*lelv*ldim)

      common /scrvh/ bmsk(lx1*ly1*lz1*lelv)
     $             , u1  (lx1*ly1*lz1*lelv)

      common /scrmg/ r1  (lx1*ly1*lz1*lelv)
     $             , r2  (lx1*ly1*lz1*lelv)
     $             , r3  (lx1*ly1*lz1*lelv)
     $             , r4  (lx1*ly1*lz1*lelv)
      
c
      nelc = nelv            ! number of elements in convecting field
      if (ifield.eq.ifldmhd) nelc = nelfld(ifield)

      nc  = cs(0)            ! number of stored convecting fields

      ln  = lx1*ly1*lz1*lelt
      n   = nx1*ny1*nz1*nelfld(ifield)
      m   = nxd*nyd*nzd*nelc*ndim

      if (ifield.eq.ifldmhd) then
         call col3(bmsk,bintm1,msk,n)
      elseif (ifield.eq.1) then
         call col3(bmsk,binvm1,msk,n)
      else ! if (ifield.eq.2) then
         call col3(bmsk,bintm1,msk,n)
      endif

      call char_conv1
     $   (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs(1),nc,ct,u1,r1,r2,r3,r4)

      return
      end
c-----------------------------------------------------------------------
      subroutine char_conv1
     $   (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs,nc,ct,u1,r1,r2,r3,r4)

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'

      real    p0(n),u(n),ulag(ln,1),bmsk(n),c(m,0:nc),cs(0:nc)

      real    ct(m),u1(n),r1(n),r2(n),r3(n),r4(n) ! work arrays

      integer gsl


!     Convect over last NBD steps using characteristics scheme

!              n-q                                      n-1
!     Given u(t    , X ),  q = 1,2,...,nbd, compute  phi     ( := p0 )

!        n-1       nbd   ~n-q
!     phi     :=  sum    u
!                  q=1

!          ~n-q             ~n-q
!     each u     satisfies  u    := v  such that


!     dv
!     -- + C.grad v = 0  t \in [t^n-q,t^n],   v(t^n-q,X) = u(t^n-q,X)
!     dt


      tau = time-vlsum(dtlag,nbd)         ! initialize time for u^n-k
      call int_vel (ct,tau,c,m,nc,cs,nid) ! ct(t) = sum w_k c(.,k)

      call rzero(p0,n)

      do ilag = nbd,1,-1

         um = 0
         if (ilag.eq.1) then
            do i=1,n
               p0(i) = p0(i)+bd(ilag+1)*u(i)
               um=max(um,u(i))
            enddo
         else
            do i=1,n
               p0(i) = p0(i)+bd(ilag+1)*ulag(i,ilag-1)
               um=max(um,ulag(i,ilag-1))
            enddo
         endif

c        write(6,1) istep,ilag,bd(ilag),bd(ilag+1),um
c 1      format(i5,i4,1p3e14.5,' bdf')

         dtau = dtlag(ilag)/ntaubd
         do itau = 1,ntaubd ! ntaubd=number of RK4 substeps (typ. 1 or 2)

            tau1 = tau + dtau

            c1 = 1.
            c2 = -dtau/2.
            c3 = -dtau
            th = tau+dtau/2.

            call conv_rhs(r1,p0,ct,bmsk,gsl)         !  STAGE 1

            call int_vel (ct,th,c,m,nc,cs,nid)       !  STAGE 2
            call add3s12 (u1,p0,r1,c1,c2,n)
            call conv_rhs(r2,u1,ct,bmsk,gsl)

            call add3s12 (u1,p0,r2,c1,c2,n)          !  STAGE 3
            call conv_rhs(r3,u1,ct,bmsk,gsl)

            call int_vel (ct,tau1,c,m,nc,cs,nid)     !  STAGE 4
            call add3s12 (u1,p0,r3,c1,c3,n)
            call conv_rhs(r4,u1,ct,bmsk,gsl)

            c1 = -dtau/6.
            c2 = -dtau/3.
            do i=1,n
               p0(i) = p0(i)+c1*(r1(i)+r4(i))+c2*(r2(i)+r3(i))
            enddo
            tau = tau1
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine int_vel(c_t,t0,c,n,nc,ct,nid)

c     Interpolate convecting velocity field c(t_1,...,t_nconv) to
c     time t0 and return result in c_t.

c     Ouput:   c_t = sum wt_k * ct_i(k) 

c     Here, t0 is the time of interest

      real c_t(n),c(n,0:nc),ct(0:nc)
c
      parameter (lwtmax=10)
      real wt(0:lwtmax)
c
      if (nc.gt.lwtmax) then
         write(6,*) nid,'ERROR int_vel: lwtmax too small',lwtmax,m0
         call exitt
      endif
c
      no = nc-1
      call fd_weights_full(t0,ct(0),no,0,wt)  ! interpolation weights

      call rzero(c_t,n)
      do j=1,n
      do i=0,no
         c_t(j) = c_t(j) + wt(i)*c(j,i)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine conv_rhs (du,u,c,bmsk,gsl)
c
      include 'SIZE'
      include 'TOTAL'
c
c     apply convecting field c(1,ndim) to scalar field u(1)
c
      real du(1),u(1),c(1),bmsk(1)
      integer gsl
c
      logical ifconv
c
c     ifconv = .false.
      ifconv = .true.
c
      n = nx1*ny1*nz1*nelv
c
      if (ifconv) then

         if (ifcons) then
           if (if3d     ) call convop_cons_3d (du,u,c,nx1,nxd,nelv)
           if (.not.if3d) call convop_cons_2d (du,u,c,nx1,nxd,nelv)
         else
           if (if3d     ) call convop_fst_3d  (du,u,c,nx1,nxd,nelv)
           if (.not.if3d) call convop_fst_2d  (du,u,c,nx1,nxd,nelv)
         endif

         call gs_op(gsl,du,1,1,0)  !  +

      else
         call rzero   (du,n)
         return
      endif
c
      do i=1,n
         du(i) = bmsk(i)*du(i)  ! Binv * msk
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine convop_fst_3d(du,u,c,mx,md,nel)
c
      include 'SIZE'
c
c     apply convecting field c to scalar field u
c
      real du(mx*mx*mx,nel)
      real  u(mx*mx*mx,nel)
      real  c(md*md*md,nel,3)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd)
c
      logical if3d,ifd
      integer e
c
      if3d = .true.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      ifd=.true.
c
      nrstd = md**3
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convop_fst')
c
      do e=1,nel
         call grad_rstd(ur,us,ut,u(1,e),mx,md,if3d,ud)
         if (ifd) then    ! dealiased
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
               ud(i) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)+c(i,e,3)*ut(i)
            enddo
            call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
         else
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
               du(i,e) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)+c(i,e,3)*ut(i)
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine convop_fst_2d(du,u,c,mx,md,nel)
c
      include 'SIZE'
c
c     apply convecting field c to scalar field u
c
      real du(mx*mx,nel)
      real  u(mx*mx,nel)
      real  c(md*md,nel,2)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd)
c
      logical if3d,ifd
      integer e
c
      if3d = .false.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      ifd=.true.
c
      nrstd = md**2
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convop_fst')
c
      do e=1,nel
         call grad_rstd(ur,us,ut,u(1,e),mx,md,if3d,ud)
         if (ifd) then    ! dealiased
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
               ud(i) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)
            enddo
            call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
         else
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
               du(i,e) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rstd(ur,us,ut,u,mx,md,if3d,ju) ! GLL->GL grad

      include 'SIZE'
      include 'DXYZ'

      real    ur(1),us(1),ut(1),u(1),ju(1)
      logical if3d

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      call intp_rstd(ju,u,mx,md,if3d,0) ! 0 = forward

      m0 = md-1
      call get_dgl_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,ju,m0,1,dg(ip),dgt(ip))
      else
         call local_grad2(ur,us   ,ju,m0,1,dg(ip),dgt(ip))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd(ju,u,mx,md,if3d,idir) ! GLL->GL interpolation

c     GLL interpolation from mx to md.

c     If idir ^= 0, then apply transpose operator  (md to mx)

      include 'SIZE'

      real    ju(1),u(1)
      logical if3d

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      parameter (ld=2*lxd)
      common /ctmp0/ w(ld**ldim,2)

      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')

      ldw = 2*(ld**ldim)

      call get_int_ptr (i,mx,md)
c
      if (idir.eq.0) then
         call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else
         call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_int(jgl,jgt,mp,np,w)
c
c     Generate interpolation from np GLL points to mp GL points
c
c        jgl  = interpolation matrix, mapping from velocity nodes to pressure
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
      real jgl(mp,np),jgt(np*mp),w(1)
c
      iz = 1
      id = iz + np
c
      call zwgll (w(iz),jgt,np)
      call zwgl  (w(id),jgt,mp)
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo
c
      call transpose(jgt,np,jgl,mp)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_dgl(dgl,dgt,mp,np,w)
c
c     Generate derivative from np GL points onto mp GL points
c
c        dgl  = interpolation matrix, mapping from velocity nodes to pressure
c        dgt  = transpose of interpolation matrix
c        w    = work array of size (3*np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
c
      real dgl(mp,np),dgt(np*mp),w(1)
c
c
      iz = 1
      id = iz + np
c
      call zwgl  (w(iz),dgt,np)  ! GL points
      call zwgl  (w(id),dgt,mp)  ! GL points
c
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
         enddo
      enddo
c
      call transpose(dgt,np,dgl,mp)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
      include 'SIZE'            ! need nid
      character*5  avar5,lvar5
      character*10 sub_name10

      if (n.gt.m) then
         write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i12,2(1x,a5),1x,a10)
         call exitti('lim_chk problem. $',n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md) ! GLL-->GL pointer

c     Get pointer to jgl() for interpolation pair (mx,md)

      include 'SIZE'

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c
      if (ip.eq.0) then
c
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
c
         call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgl_ptr (ip,mx,md)
c
c     Get pointer to GL-GL interpolation dgl() for pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c
      parameter (ld=2*lxd)
      common /jgrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pdg (ij)

      if (ip.eq.0) then

         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c
         call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_char(ct,c,ux,uy,uz,nelc,tau,ifnew)
      include 'SIZE'
      include 'TSTEP'

      real ct(0:1)               ! time stamps for saved field (0=#flds)
      real c(1)                  ! saved vel. fields, dealiased etc.
      real ux(1),uy(1),uz(1)     ! input vel. field
      integer nelc               ! number of elements in conv. field
      logical ifnew              ! =true if shifting stack of fields

      numr      = lxd*lyd*lzd*lelv*ldim*(lorder+1)
      denr      = nxd*nyd*nzd*nelv*ndim
      nconv_max = numr/denr
      if (nconv_max.lt.nbdinp+1) 
     $   call exitti(
     $     'ABORT: not enough memory for characteristics scheme!$',
     $     nconv_max)

      nc = ct(0)

      m  = nxd*nyd*nzd*nelc*ndim

c     write(6,*) nelc,ifnew,' set conv_char',istep,nc,nconv_max
      call set_ct_cvx
     $    (ct,c,m,ux,uy,uz,tau,nc,nconv_max,nelc,ifnew)

      nc = min (nc,nbdinp)
      ct(0) = nc  ! store current count

      return
      end
c-----------------------------------------------------------------------
      subroutine set_ct_cvx(ct,c,m,u,v,w,tau,nc,mc,nelc,ifnew)
      include 'SIZE'
      include 'INPUT'  ! ifcons

      real ct(0:1),c(m,1)
      real u(1),v(1),w(1)
      logical ifnew

      if (ifnew) then

c        Shift existing convecting fields
c        Note:  "1" entry is most recent

         nc = nc+1
         nc = min(nc,mc)
         ct(0) = nc

         do i=nc,2,-1
            call copy(c(1,i),c(1,i-1),m)
            ct(i) = ct(i-1)
         enddo
      endif

c     Save time and map the current velocity to rst coordinates.

      ix = 1
      iy = ix + nxd*nyd*nzd*nelc
      iz = iy + nxd*nyd*nzd*nelc

      if (ifcons) then
         call set_convect_cons(c(ix,1),c(iy,1),c(iz,1),u,v,w)
      else
         call set_convect_new (c(ix,1),c(iy,1),c(iz,1),u,v,w)
      endif

      ct(1) = tau

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rst(ur,us,ut,u,md,if3d) ! Gauss-->Gauss grad

      include 'SIZE'
      include 'DXYZ'

      real    ur(1),us(1),ut(1),u(1)
      logical if3d

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      m0 = md-1
      call get_dgl_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,u,m0,1,dg(ip),dgt(ip))
      else
         call local_grad2(ur,us   ,u,m0,1,dg(ip),dgt(ip))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine convect_new(bdu,u,ifuf,cx,cy,cz,ifcf)

C     Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
C
      include 'SIZE'
      include 'TOTAL'

      real bdu(1),u(1),cx(1),cy(1),cz(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e

      call set_dealias_rx

      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd

      nxyzu = nxyz1
      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
      if (ifcf) nxyzc = nxyzd

      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ib = 1    ! pointer to scalar field Bdu


      do e=1,nelv

         if (ifcf) then

            call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
            call copy(tr(1,2),cy(ic),nxyzd)
            if (if3d) call copy(tr(1,3),cz(ic),nxyzd)

         else  ! map coarse velocity to fine mesh (C-->F)

           call intp_rstd(fx,cx(ic),nx1,nxd,if3d,0) ! 0 --> forward
           call intp_rstd(fy,cy(ic),nx1,nxd,if3d,0) ! 0 --> forward
           if (if3d) call intp_rstd(fz,cz(ic),nx1,nxd,if3d,0) ! 0 --> forward

           if (if3d) then  ! Convert convector F to r-s-t coordinates

             do i=1,nxyzd
               tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
               tr(i,2)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
               tr(i,3)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
             enddo

           else

             do i=1,nxyzd
               tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
               tr(i,2)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
             enddo

           endif

         endif

         if (ifuf) then
            call grad_rst(ur,us,ut,u(iu),nxd,if3d)
         else
            call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward
            call grad_rst(ur,us,ut,uf,nxd,if3d)
         endif

         if (if3d) then
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
            enddo
         else
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)
            enddo
         endif
         call intp_rstd(bdu(ib),uf,nx1,nxd,if3d,1) ! Project back to coarse

         ic = ic + nxyzc
         iu = iu + nxyzu
         ib = ib + nxyz1

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine convect_cons(bdu,u,ifuf,cx,cy,cz,ifcf)

c     Compute dealiased form:  J^T Bf *div. JC Ju w/ correct Jacobians

c     conservative form


      include 'SIZE'
      include 'TOTAL'

      real bdu(1),u(1),cx(1),cy(1),cz(1)

      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ uf(ltd),cf(ltd),cu(ltd)
     $             , cr(ltd),cs(ltd),ct(ltd)


      integer e

      call set_dealias_rx

      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd

      nxyzu = nxyz1
      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
      if (ifcf) nxyzc = nxyzd

      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ib = 1    ! pointer to scalar field Bdu

      do e=1,nelv

        call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward

        call rzero(cu,nxyzd)
        do i=1,ndim

         if (ifcf) then  ! C is already on fine mesh

           call exitt  ! exit for now

         else  ! map coarse velocity to fine mesh (C-->F)

           if (i.eq.1) call intp_rstd(cf,cx(ic),nx1,nxd,if3d,0) ! 0 --> forward
           if (i.eq.2) call intp_rstd(cf,cy(ic),nx1,nxd,if3d,0) ! 0 --> forward
           if (i.eq.3) call intp_rstd(cf,cz(ic),nx1,nxd,if3d,0) ! 0 --> forward

           call col2(cf,uf,nxyzd)   !  collocate C and u on fine mesh

           call grad_rst(cr,cs,ct,cf,nxd,if3d)  ! d/dr (C_i*u)

           if (if3d) then

             do j=1,nxyzd
               cu(j)=cu(j)
     $              +cr(j)*rx(j,i,e)+cs(j)*rx(j,i+3,e)+ct(j)*rx(j,i+6,e)
             enddo

           else  ! 2D

             do j=1,nxyzd
               cu(j)=cu(j)
     $              +cr(j)*rx(j,i,e)+cs(j)*rx(j,i+2,e)
             enddo

           endif
         endif
        enddo

        call intp_rstd(bdu(ib),cu,nx1,nxd,if3d,1) ! Project back to coarse

        ic = ic + nxyzc
        iu = iu + nxyzu
        ib = ib + nxyz1

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_convect_cons(cx,cy,cz,ux,uy,uz)

c     Put vx,vy,vz on fine mesh (for conservation form)


      include 'SIZE'
      include 'TOTAL'

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real cx(ltd,1),cy(ltd,1),cz(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)

      integer e

      call set_dealias_rx

      do e=1,nelv    ! Map coarse velocity to fine mesh (C-->F)

         call intp_rstd(cx(1,e),ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         call intp_rstd(cy(1,e),uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(cz(1,e),uz(1,e),nx1,nxd,if3d,0) ! 0 --> forward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_convect_new(cr,cs,ct,ux,uy,uz)
C
C     Put vxd,vyd,vzd into rst form on fine mesh
C
C     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
C
      include 'SIZE'
      include 'TOTAL'

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real cr(ltd,1),cs(ltd,1),ct(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)

      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e

      call set_dealias_rx

      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd

      ic = 1    ! pointer to vector field C

      do e=1,nelv 

c        Map coarse velocity to fine mesh (C-->F)

         call intp_rstd(fx,ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         call intp_rstd(fy,uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(fz,uz(1,e),nx1,nxd,if3d,0) ! 0 --> forward

c        Convert convector F to r-s-t coordinates

         if (if3d) then

           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
              cs(i,e)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
              ct(i,e)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
           enddo

         else

           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
              cs(i,e)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
           enddo

         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_char_mask(mask,u,v,w) ! mask for hyperbolic system 

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'TSTEP'
      integer msk(0:1)
      character      cb*3
      parameter (lxyz1=lx1*ly1*lz1)
      common /ctmp1/ work(lxyz1,lelt)

      real mask(lxyz1,1),u(lxyz1,1),v(lxyz1,1),w(lxyz1,1)

      integer e,f

      nfaces= 2*ndim
      ntot1 = nx1*ny1*nz1*nelv
      call rzero (work,ntot1)
      call rone  (mask,NTOT1)

      ifldv = 1
      do 100 e=1,nelv
      do 100 f=1,nfaces
         cb=cbc(f,e,ifldv)
         if (cb(1:1).eq.'v' .or. cb(1:1).eq.'V') then

           call faccl3 (work(1,e),u(1,e),unx(1,1,f,e),f)
           call faddcl3(work(1,e),v(1,e),uny(1,1,f,e),f)
           if (if3d) 
     $     call faddcl3(work(1,e),w(1,e),unz(1,1,f,e),f)

           call fcaver (vaver,work,e,f)

           if (vaver.lt.0) call facev (mask,e,f,0.0,nx1,ny1,nz1)
         endif
         if (cb(1:2).eq.'ws' .or. cb(1:2).eq.'WS') 
     $   call facev (mask,e,f,0.0,nx1,ny1,nz1)
 100  continue

      return
      end
c-----------------------------------------------------------------------
      subroutine advchar
c
c     Compute convective contribution using 
c     operator-integrator-factor method (characteristics).
c
      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'PARALLEL'

      include 'CTIMER'

      common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()

      common /scruz/ phx  (lx1*ly1*lz1*lelt)
     $ ,             phy  (lx1*ly1*lz1*lelt)
     $ ,             phz  (lx1*ly1*lz1*lelt)
     $ ,             hmsk (lx1*ly1*lz1*lelt)

#ifndef NOTIMER
      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()
#endif

      dti = 1./dt
      n   = nx1*ny1*nz1*nelv

      call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system 

      call char_conv(phx,vx,vxlag,hmsk,c_vx,ct_vx,gsh_fld(1))
      call char_conv(phy,vy,vylag,hmsk,c_vx,ct_vx,gsh_fld(1))

      if (if3d) then

        call char_conv(phz,vz,vzlag,hmsk,c_vx,ct_vx,gsh_fld(1))

        do i=1,n
           h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
           bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
           bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
           bfz(i,1,1,1) = bfz(i,1,1,1)+phz(i)*h2i
        enddo

      else

        do i=1,n
           h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
           bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
           bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
        enddo

      endif

#ifndef NOTIMER
      tadvc=tadvc+(dnekclock()-etime1)
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine convch

c     Compute convective contribution using 
c     operator-integrator-factor method (characteristics).

      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'PARALLEL'
      include 'CTIMER'

      common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()

      common /scruz/ phi  (lx1*ly1*lz1*lelt)
     $ ,             hmsk (lx1*ly1*lz1*lelt)

      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()

      n   = nx1*ny1*nz1*nelv
      dti = 1./dt

      if (ifield.eq.2) then  ! set convecting velocity and mask
c        call setup_convect(1)
         call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system 
      endif


      call char_conv(phi,t(1,1,1,1,ifield-1),tlag(1,1,1,1,1,ifield-1)
     $        ,hmsk,c_vx,ct_vx,gsh_fld(1))

c     pmax = glamax(phi,n)
c     qmax = glamax(vtrans(1,1,1,1,2),n)
c     write(6,*) istep,dti,pmax,' pmax'

      do i=1,n
         bq(i,1,1,1,ifield-1) = bq(i,1,1,1,ifield-1)
     $          + phi(i)*bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)*dti
      enddo

      tadvc=tadvc+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------
      subroutine convop_cons_3d(du,u,c,mx,md,nel) ! Conservation form

c     Apply convecting field c to scalar field u, conservation form d/dxj cj phi

c     Assumes that current convecting field is on dealias mesh, in c()

      include 'SIZE'
      include 'DEALIAS'
      include 'GEOM'

      real du(mx*mx*mx,nel)
      real  u(mx*mx*mx,nel)
      real  c(md*md*md,nel,3)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju

      logical if3d,ifd
      integer e

      if3d = .true.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      ifd=.true.

      nrstd = md**3
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convp_cons')

      do e=1,nel

         call intp_rstd (ju,u(1,e),mx,md,if3d,0) ! 0 = forward; on Gauss points!
         call rzero     (ud,nrstd)

         do j=1,ndim
            do i=1,nrstd
               tu(i)=c(i,e,j)*ju(i)   ! C_j*T
            enddo
            call grad_rst(ur,us,ut,tu,md,if3d)  ! Already on fine (Gauss) mesh

            j0 = j+0
            j3 = j+3
            j6 = j+6
            do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
               ud(i)=ud(i)
     $              +rx(i,j0,e)*ur(i)+rx(i,j3,e)*us(i)+rx(i,j6,e)*ut(i)
            enddo
         enddo

         call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine convop_cons_2d(du,u,c,mx,md,nel) ! Conservation form

c     Apply convecting field c to scalar field u, conservation form d/dxj cj phi

c     Assumes that current convecting field is on dealias mesh, in c()

      include 'SIZE'
      include 'GEOM'
      include 'TSTEP'


      real du(mx*mx,nel)
      real  u(mx*mx,nel)
      real  c(md*md,nel,2)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju

      logical if3d,ifd
      integer e

      if3d = .false.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.

      nrstd = md**2
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convp_cons')

      if (nio.eq.0.and.istep.lt.3) write(6,*) 'convp_cons',istep

      do e=1,nel

         call intp_rstd (ju,u(1,e),mx,md,if3d,0) ! 0 = forward; on Gauss points!
         call rzero     (ud,nrstd)

c        call outmat(c(1,e,1),md,md,'fine u',e)
c        call outmat(c(1,e,2),md,md,'fine v',e)
c        call outmat(ju      ,md,md,'fine T',e)

         do j=1,ndim
            do i=1,nrstd
               tu(i)=c(i,e,j)*ju(i)   ! C_j*T
            enddo
            call grad_rst(ur,us,ut,tu,md,if3d)  ! Already on fine (Gauss) mesh

            j0 = j+0
            j2 = j+2
            do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
               ud(i)=ud(i)+rx(i,j0,e)*ur(i)+rx(i,j2,e)*us(i)
            enddo
         enddo

         call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose

      enddo

c     call exitti('convop_cons_2d$',istep)

      return
      end
c-----------------------------------------------------------------------
