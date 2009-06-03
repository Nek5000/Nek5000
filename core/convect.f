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
      nc  = cs(0)                 ! Number of stored convecting fields

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




      tau = time-vlsum(dtlag,nbd) ! initialize time for u^n-k

      call int_vel (ct,tau,c,m,nc,cs,nid) ! ct(t) = sum w_k c(.,k)
      call rzero(p0,n)
c
      do ilag = nbd,1,-1

         if (ilag.eq.1) then
            do i=1,n
               p0(i) = p0(i)+bd(ilag+1)*u(i)
            enddo
         else
            do i=1,n
               p0(i) = p0(i)+bd(ilag+1)*ulag(i,ilag-1)
            enddo
         endif

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
c
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
c
      call rzero(c_t,n)
      k = 0
      do i=0,no
         do j=1,n
            c_t(j) = c_t(j) + wt(k)*c(j,i)
         enddo
         k = k+1
      enddo
c
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
         if (if3d) then
            call convop_fst_3d (du,u,c,nx1,nxd,nelv)
         else
            call convop_fst_2d (du,u,c,nx1,nxd,nelv)
         endif
         call gs_op(gsl,du,1,1,0)
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
      subroutine grad_rstd(ur,us,ut,u,mx,md,if3d,ju)
c
      include 'SIZE'
      include 'DXYZ'
c
      real    ur(1),us(1),ut(1),u(1),ju(1)
      logical if3d
c
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      m0 = md-1
      call intp_rstd(ju,u,mx,md,if3d,0) ! 0 = forward

      call get_dgl_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,ju,m0,1,dg(ip),dgt(ip))
      else
         call local_grad2(ur,us   ,ju,m0,1,dg(ip),dgt(ip))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      logical if3d
c
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c
      parameter (ld=2*lxd)
      common /ctmp0/ w(ld**ldim,2)
c
      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
c
      ldw = 2*(ld**ldim)
c
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
c
      if (n.gt.m) then
         write(6,*) nid,n,m
c        write(6,*) nid,n,m,avar5,lvar5,sub_name10
c        write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i9,2(1x,a5),1x,a10)
         call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md)
c
c     Get pointer to jgl() for interpolation pair (mx,md)
c
      include 'SIZE'
c
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
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pdg (ij)
c
      if (ip.eq.0) then
c
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
      subroutine set_conv(c,ux,uy,uz,nelc,ct,tau)
      include 'SIZE'
      include 'TSTEP'
c
      real c(1)
      real ct(0:1)
      real ux(1),uy(1),uz(1)

      numr      = lxd*lyd*lzd*lelv*ldim*(lorder-1)
      denr      = nxd*nyd*nzd*nelv*ndim
      nconv_max = numr/denr
      nconv_max = min(nconv_max,nbdinp+1)

      nc = ct(0)
      m  = nxd*nyd*nzd*nelc*ndim

      call set_conv_tau
     $    (c,m,ux,uy,uz,ct(1),tau,nc,nconv_max,nx1,nxd,nelc)

      nc = min (nc,nbdinp)
      ct(0) = nc  ! store current count
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_tau(c,m,u,v,w,ct,tau,nc,mc,mx,md,nelc)
c
      real c(m,1),ct(1)
      real u(1),v(1),w(1)
c
c     First, shift existing convecting fields
c
c     Note:  "1" entry is most recent
c
      nc = nc+1
      nc = min(nc,mc)
c
      do i=nc,2,-1
         call copy(c(1,i),c(1,i-1),m)
         ct(i) = ct(i-1)
      enddo
c
c     Next, save time and map the current velocity to rst coordinates.
c
      call set_conv_fld(c(1,1),u,v,w,md,mx,nelc)
      ct(1) = tau
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_fld(c,u,v,w,md,mx,nel)
c
c     Set up convecting field C
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
c
      real c(1),u(1),v(1),w(1)
c
      if (mx.eq.md) then
         n = nel*mx**ndim
         call set_conv_fld_alias(c,u,v,w,n)
      else
         if (if3d) then
            call set_conv_fld_dealias_3d(c,u,v,w,md,mx,nel)
         else
            call set_conv_fld_dealias_2d(c,u,v,md,mx,nel)
         endif
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_fld_alias(c,u,v,w,n)
c
c     Set up convecting field c, no dealias
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
c
      real c(n,3),u(n),v(n),w(n)
c
      write(6,*) '... rx * mass matrix ???  ERROR! '
      call exitt
      if (if3d) then
         do i=1,n
            c(i,1) = ( rxm1(i,1,1,1)*u(i) + rym1(i,1,1,1)*v(i)
     $             +   rzm1(i,1,1,1)*w(i))* bm1 (i,1,1,1)
            c(i,2) = ( sxm1(i,1,1,1)*u(i) + sym1(i,1,1,1)*v(i)
     $             +   szm1(i,1,1,1)*w(i))* bm1 (i,1,1,1)
            c(i,3) = ( txm1(i,1,1,1)*u(i) + tym1(i,1,1,1)*v(i)
     $             +   tzm1(i,1,1,1)*w(i))* bm1 (i,1,1,1)
         enddo
      elseif (ifaxis) then
         write(6,*) 'no axis in set_conv_fld'
         call exitt
      else
         do i=1,n
            c(i,1) = (rxm1(i,1,1,1)*u(i)+rym1(i,1,1,1)*v(i))
     $             *   bm1(i,1,1,1)      !     *0  !DEBUG
            c(i,2) = (sxm1(i,1,1,1)*u(i)+sym1(i,1,1,1)*v(i))
     $             *   bm1(i,1,1,1)      !     *0  !DEBUG
         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_fld_dealias_3d(c,u,v,w,md,mx,nel)
c
c     Set up convecting field c, dealias
c
c     Currently, we're using the poor man's version of dealiasing,
c     where we cheat on the Jacobian and metric term.   That is,
c     our variational form reads:
c
c                   T
c     (v,Cu) := (Jv)  B J(C . X ) D  Ju                    (1)
c                              r   r
c               
c     as opposed to:
c
c                   T
c     (v,Cu) := (Jv)  B JC . JX D Ju                       (2)
c                              r r
c               
c
c     In form (1), which is coded below, we are collocating the convecting 
c     field C with the Jacobian-metric product, Xr, prior to interpolating 
c     onto the dealiasing mesh.   Form (2) is more correct, but would require
c     either storage or recompuation of the interpolated tensor JXr.
c
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
c
      real c(md*md*md,nel,3)
      real u(mx*mx*mx,nel),v(mx*mx*mx,nel),w(mx*mx*mx,nel)
c
      integer mo,e
      save    mo
      data    mo / 0 /
c
      common /wdealias/ wgld(lxd),wwwd(lxd*lxd*lxd)
c
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ cc(ldd)
c
      if (md.ne.mo) then
         call zwgl(c,wgld,md)
         mo = md
c
         l = 0
         do k=1,md
         do j=1,md
         do i=1,md
            l=l+1
            wwwd(l) = wgld(i)*wgld(j)*wgld(k)
         enddo
         enddo
         enddo
      endif
c
      nxyd = md**ndim
      nxyz = mx**ndim
      do e=1,nel
c
         do i=1,nxyz
            cc(i) = rxm1(i,1,1,e)*u(i,e) + rym1(i,1,1,e)*v(i,e)
     $            + rzm1(i,1,1,e)*w(i,e)
         enddo
         call intp_rstd(c(1,e,1),cc,mx,md,if3d,0) ! 0 --> forward
         do i=1,nxyd
            c(i,e,1) = wwwd(i)*c(i,e,1)
         enddo
c
         do i=1,nxyz
            cc(i) = sxm1(i,1,1,e)*u(i,e) + sym1(i,1,1,e)*v(i,e)
     $            + szm1(i,1,1,e)*w(i,e)
         enddo
         call intp_rstd(c(1,e,2),cc,mx,md,if3d,0) ! 0 --> forward
         do i=1,nxyd
            c(i,e,2) = wwwd(i)*c(i,e,2)
         enddo
c
         do i=1,nxyz
            cc(i) = txm1(i,1,1,e)*u(i,e) + tym1(i,1,1,e)*v(i,e)
     $            + tzm1(i,1,1,e)*w(i,e)
         enddo
         call intp_rstd(c(1,e,3),cc,mx,md,if3d,0) ! 0 --> forward
         do i=1,nxyd
            c(i,e,3) = wwwd(i)*c(i,e,3)
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_fld_dealias_2d(c,u,v,md,mx,nel)
c
c     Set up convecting field c, dealias
c
c     Currently, we're using the poor man's version of dealiasing,
c     where we cheat on the Jacobian and metric term.   That is,
c     our variational form reads:
c
c                   T
c     (v,Cu) := (Jv)  B J(C . X ) D  Ju                    (1)
c                              r   r
c               
c     as opposed to:
c
c                   T
c     (v,Cu) := (Jv)  B JC . JX D Ju                       (2)
c                              r r
c               
c
c     In form (1), which is coded below, we are collocating the convecting 
c     field C with the Jacobian-metric product, Xr, prior to interpolating 
c     onto the dealiasing mesh.   Form (2) is more correct, but would require
c     either storage or recompuation of the interpolated tensor JXr.
c
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
c
      real c(md*md,nel,2)
      real u(mx*mx,nel),v(mx*mx,nel)
c
      integer mo,e
      save    mo
      data    mo / 0 /
c
      common /wdealias/ wgld(lxd),wwwd(lxd*lxd*lxd)
c
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ cc(ldd)
c
      if (md.ne.mo) then
         call zwgl(c,wgld,md)
         mo = md
c
         l = 0
         do j=1,md
         do i=1,md
            l=l+1
            wwwd(l) = wgld(i)*wgld(j)
         enddo
         enddo
      endif
c
      nxyd = md**ndim
      nxyz = mx**ndim
      do e=1,nel
c
         do i=1,nxyz
            cc(i) = rxm1(i,1,1,e)*u(i,e) + rym1(i,1,1,e)*v(i,e)
         enddo
         call intp_rstd(c(1,e,1),cc,mx,md,if3d,0) ! 0 --> forward
         do i=1,nxyd
            c(i,e,1) = wwwd(i)*c(i,e,1)
         enddo
c
         do i=1,nxyz
            cc(i) = sxm1(i,1,1,e)*u(i,e) + sym1(i,1,1,e)*v(i,e)
         enddo
         call intp_rstd(c(1,e,2),cc,mx,md,if3d,0) ! 0 --> forward
         do i=1,nxyd
            c(i,e,2) = wwwd(i)*c(i,e,2)
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rst(ur,us,ut,u,md,if3d)
c
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
C
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
      subroutine set_convect_new(cr,cs,ct)
C
C     Put vxd,vyd,vzd into rst form on fine mesh
C
C     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
C
      include 'SIZE'
      include 'TOTAL'

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real cr(ltd,1),cs(ltd,1),ct(ltd,1)

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

         call intp_rstd(fx,vx(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> forward
         call intp_rstd(fy,vy(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(fz,vz(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> forward

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
