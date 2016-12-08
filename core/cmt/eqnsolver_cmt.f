      subroutine assemble_h(e,eq_num)
      include 'SIZE'
      include 'CMTDATA'
      include 'INPUT'

      integer  e,eq_num

!     !This subroutine will compute the convective and diffusive 
!     !components of h
      call rzero(totalh,3*lxd*lyd*lzd)
      if (nxd.gt.nx1) then
         call evaluate_dealiased_conv_h(e,eq_num)
      else
         call evaluate_aliased_conv_h(e,eq_num)
      endif
      if (ifvisc.and.eq_num .gt. 1) call evaluate_diff_h(e,eq_num)
      call add_conv_diff_h
      
      return
      end

!-----------------------------------------------------------------------

      subroutine evaluate_dealiased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'
     
      integer  e,eq

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ju1(ldd),ju2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju1,ju2

      n=nxd*nyd*nzd

c     if (eq .eq. 1) then ! convective flux of mass=rho u_j=U_{j+1}

c        do j=1,ndim
c           call intp_rstd(convh(1,j),u(1,1,1,eq+j,e),nx1,nxd,if3d,0)
c        enddo

c     else
c To be consistent with momentum equation, for mass balance flux vector is 
c computed by multiplying rho by u_j
         call intp_rstd(ju1,phig(1,1,1,e),nx1,nxd,if3d,0)
         call intp_rstd(ju2,pr(1,1,1,e),nx1,nxd,if3d,0)

         if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j

            call intp_rstd(convh(1,1),u(1,1,1,eq,e),nx1,nxd,if3d,0)
            do j=2,ndim
               call copy(convh(1,j),convh(1,1),n)
            enddo
            call col2(convh(1,1),vxd(1,1,1,e),n)
            call col2(convh(1,2),vyd(1,1,1,e),n)
            if (if3d) call col2(convh(1,3),vzd(1,1,1,e),n)
            call add2col2(convh(1,eq-1),ju1,ju2,n)

         elseif (eq .eq. 5) then

            call intp_rstd(convh(1,1),u(1,1,1,eq,e),nx1,nxd,if3d,0)
            call add2col2(convh(1,1),ju1,ju2,n)
            do j=2,ndim
               call copy(convh(1,j),convh(1,1),n)
            enddo
            call col2(convh(1,1),vxd(1,1,1,e),n)
            call col2(convh(1,2),vyd(1,1,1,e),n)
            call col2(convh(1,3),vzd(1,1,1,e),n)

         else
            if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
            if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
            call exitt
         endif

c     endif
     
      return
      end

!-----------------------------------------------------------------------

      subroutine evaluate_aliased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ju1(ldd),ju2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju1,ju2
      integer  e,eq

      n=nxd*nyd*nzd

      call copy(ju1,phig(1,1,1,e),n)
      call copy(ju2,pr(1,1,1,e),n)

      if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j

         call copy(convh(1,1),u(1,1,1,eq,e),n)
         do j=2,ndim
            call copy(convh(1,j),convh(1,1),n)
         enddo
         call col2(convh(1,1),vxd(1,1,1,e),n)
         call col2(convh(1,2),vyd(1,1,1,e),n)
         if (if3d) call col2(convh(1,3),vzd(1,1,1,e),n)
         if(eq. gt. 1) call add2col2(convh(1,eq-1),ju1,ju2,n)

      elseif (eq .eq. 5) then

         call copy(convh(1,1),u(1,1,1,eq,e),n)
         call add2col2(convh(1,1),ju1,ju2,n)
         do j=2,ndim
            call copy(convh(1,j),convh(1,1),n)
         enddo
         call col2(convh(1,1),vxd(1,1,1,e),n)
         call col2(convh(1,2),vyd(1,1,1,e),n)
         call col2(convh(1,3),vzd(1,1,1,e),n)

      else
         if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
         if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
         call exitt
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine evaluate_diff_h(e,eq)
      include  'SIZE'
      include  'CMTDATA'
      include  'SOLN'
      include  'DG'
      include  'INPUT'
! diagnostic
      include 'GEOM'
      include 'TSTEP'
! diagnostic

      parameter (ldd=lx1*ly1*lz1)
! JH090815 need to get rid of this
      common /primd/ pvart(ldd,lelcmt,7) ! rho,u,v,w,T,P,phi

! JH020215 Congratulations. You can't use this common block once toteq>5
      common /ctmp1/ jdu(ldd,toteq),visco(lx1,ly1,lz1)

      real jdu,visco
      integer lfq,lfc,hcsize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelcmt,
     >                   hcsize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*ldim*lfq)
      common /CMTSURFLX/ oldface(hcsize),fatface(hdsize)
      integer e,eq
      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
!     real yrface(nxd*nzd,toteq) ! eh. just chuck yrface in jdu
! JH012915 I omitted diffh to save some space. totalh should be zeroed
!         before this routine is called. assemble_h will +=totalh
!         with convh.

      nxyz=nx1*ny1*nz1
      nface = 2*ndim
      nxz   = nx1*nz1
      nfq=nx1*nz1*2*ndim*nelt
      nstate=nqq

!-----------------------------------------------------------------------
! JH090815 must get rid of this
      if (e .eq. 1) then ! do the worst thing ever
         call copy(pvart(1,1,1),vtrans,nxyz*nelt)
         call copy(pvart(1,1,2),vx,nxyz*nelt)
         call copy(pvart(1,1,3),vy,nxyz*nelt)
         call copy(pvart(1,1,4),vz,nxyz*nelt)
         call copy(pvart(1,1,5),t,nxyz*nelt)
         call copy(pvart(1,1,6),pr,nxyz*nelt)
         call copy(pvart(1,1,7),phig,nxyz*nelt)
      endif
! JH090815 must get rid of that
!-----------------------------------------------------------------------

! MS050615 gradu contains the the surface integral contribution and
!          the solution to the DG formulation of Sij=gradU.
!  AuxFlux call was removed from here and placed outside the equation
!    loop, rigth after compute_gradients is called.
      do j=1,ndim
         do k=1,ndim
            ieijk=0
            if (eq .lt. 5) ieijk=eijk3(eq-1,j,k) ! does this work in 2D?

            if (ieijk .eq. 0) then
              do l=1,toteq
               call copy(jdu(1,l),gradu(1,1,1,l,k),nxyz)
              enddo
              call tauij_gdu(totalh(1,j),visco,pvart,pvart,7,jdu,jdu,
     >                 toteq,nxyz,.true.,eq,j,k,e,
     >                 vdiff(1,1,1,e,imu),
     >                 vdiff(1,1,1,e,ilam),
     >                 vdiff(1,1,1,e,iknd),
     >                 vtrans(1,1,1,e,icv))
            endif
         enddo

! strip faces for numerical flux of Hd. BR1, Sij=gradU
         call store_gdu_hstarsij(fatface,totalh(1,j),e,eq,j)
! MS050615 - if the auxiliary variable is GgradU, then the above
!            store_gdu_hstarsij call is NOT done for BR1!!!
!            It would be done ONLY for SIP (Hartmann & Houston 2008).
      enddo

! MS050615 - if the auxiliary variable is GgradU, then AuxFlux must
!            be called here unless you want gradu to be toteq times
!            bigger!
!      do j=1,ndim
!! get -area/bm1*Gijkl[[Ul]]k refined into fine faces, jdu(:,1)
!         call AuxFlux(jdu(1,1),oldface(iqm),oldface(iuj),jdu(1,2),
!     >                jdu(1,3),e,eq,j,nstate)
!! and add_face2full to turn totalh into Sij, i=eq
!         call add_face2full_cmt(1,nxd,nyd,nzd,iface_fine(1,e),
!     >                          totalh(1,j),jdu(1,1))
!         if (ifbr1)
!     >   call store_gdu_hstarsij(fatface,totalh(1,j),e,eq,j)
!      enddo

      call cmult(totalh,-1.0,3*nxyz)

      return
      end

!-----------------------------------------------------------------------

      subroutine add_conv_diff_h
      include  'SIZE'
      include  'CMTDATA'

      n = nxd*nyd*nzd*3

!BLAS is also a possibility here

!     call add3(totalh,convh,diffh,n)
      call add2(totalh,convh,n)

      return
      end

!----------------------

      subroutine flux_div_integral_dealiased(e,eq)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      integer  dir
      parameter (ldd=lxd*lyd*lzd)
      parameter (ldg=lxd**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      nrstd=ldd
      nxyz=nx1*ny1*nz1
      call get_dgl_ptr(ip,nxd,nxd) ! fills dg, dgt
      mdm1=nxd-1

      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      j0=0
      do j=1,ndim
         j0=j0+1
         call add2col2(ur,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      do j=1,ndim
         j0=j0+1
         call add2col2(us,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      if (if3d) then
         do j=1,ndim
            j0=j0+1
            call add2col2(ut,totalh(1,j),rx(1,j0,e),nrstd)
         enddo
         call local_grad3_t(ud,ur,us,ut,mdm1,1,dg(ip),dgt(ip),wkd)
      else
         call local_grad2_t(ud,ur,us,   mdm1,1,dg(ip),dgt(ip),wkd)
      endif

      call intp_rstd(tu,ud,nx1,nxd,if3d,1)

! multiply the residue by mass matrix. Non diagonal B should still be
! one block per element
!     call col2(ud,bm1(1,1,1,e),nxyz)  ! res = B*res  --- B=mass matrix
!     call add2(res1(1,1,1,e,eq),tu,nxyz)
! weak?
      call sub2(res1(1,1,1,e,eq),tu,nxyz)

      return
      end

!----------------------

      subroutine flux_div_integral_aliased(e,eq)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      integer  dir
      parameter (ldd=lxd*lyd*lzd)
      parameter (ldg=lxd**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      character*32 cname

      nrstd=ldd
      nxyz=nx1*ny1*nz1
      call get_dgll_ptr(ip,nxd,nxd) ! fills dg, dgt
      mdm1=nxd-1

      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      j0=0
      do j=1,ndim
         j0=j0+1
         call add2col2(ur,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      do j=1,ndim
         j0=j0+1
         call add2col2(us,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      if (if3d) then
         do j=1,ndim
            j0=j0+1
            call add2col2(ut,totalh(1,j),rx(1,j0,e),nrstd)
         enddo
         call local_grad3_t(ud,ur,us,ut,mdm1,1,d(ip),dt(ip),wkd)
      else
         call local_grad2_t(ud,ur,us,   mdm1,1,d(ip),dt(ip),wkd)
      endif

      call copy(tu,ud,nxyz)

      call sub2(res1(1,1,1,e,eq),tu,nxyz)

      return
      end

!-----------------------------------------------------------------------
      subroutine compute_forcing(e,eq_num)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'SOLN'
      include  'CMTDATA'
      include  'DEALIAS'
      
      integer e,eq_num
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      nrstd=ldd
      nxyz=nx1*ny1*nz1
      call rzero(ud,nxyz)
      if(eq_num.ne.1.and.eq_num.ne.5)then
        if(eq_num.eq.2)then
           j=1
        elseif(eq_num.eq.3)then
           j=2
        elseif(eq_num.eq.4)then
           j=2
           if(ldim.eq.3) j=3
        endif
c       write(6,*)'enter  compute_forcing ', j
        call gradl_rst(ur,us,ut,phig(1,1,1,e),lx1,if3d) ! navier1
        if (if3d) then
           j0=j+0
           j3=j+3
           j6=j+6
           do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
              ud(i)=rx(i,j0,e)*ur(i)+rx(i,j3,e)*us(i)+rx(i,j6,e)*ut(i)
           enddo
        else
           j0=j+0
           j2=j+2
           do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
              ud(i)=rx(i,j0,e)*ur(i)+rx(i,j2,e)*us(i)
           enddo
        endif
        if (eq_num.eq.4.and.ldim.eq.2)then

        else
           call col2(ud,pr(1,1,1,e),nxyz)
           call copy(convh(1,1),ud,nxyz)
           call col2(convh(1,1),jacmi(1,e),nxyz)
           call col2(convh(1,1),bm1(1,1,1,e),nxyz)  ! res = B*res
           call sub2(res1(1,1,1,e,eq_num),convh(1,1),nxyz)
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
        endif
      elseif(eq_num.eq.5)then
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cmtusrf(e)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'TSTEP'
      include 'PARALLEL'

      integer e,eg

      if(istep.eq.1)then
        n = nx1*ny1*nz1*5
        call rzero(usrf,n)
      endif
      eg = lglel(e)
      do k=1,nz1
         do j=1,ny1
            do i=1,nx1
               call NEKASGN(i,j,k,e)
               call userf(i,j,k,eg)
               usrf(i,j,k,2) = FFX
               usrf(i,j,k,3) = FFY
               usrf(i,j,k,4) = FFZ
               usrf(i,j,k,5) = (U(i,j,k,2,e)*FFX + U(i,j,k,3,e)*FFY
     &                       +  U(i,j,k,4,e)*FFZ)/ U(i,j,k,1,e)
            enddo
         enddo
      enddo

      return
      end 

!-----------------------------------------------------------------------

      subroutine penalty(qminus,qplus,ujump,e,eq,nstate)
! computes penalty function in Eq.
! of Hartmann and Houston (2008)
! This would have been about 100 times easier with the perpetual assumption
! nx1=nx1, but I decided to do some heavy lifting now.
      include 'SIZE'
      include 'CMTDATA'
      include 'GEOM'
      include 'DG'
      include 'INPUT'
      integer ctmp1pad
      parameter (ld2=lxd*lzd,ctmp1pad=(toteq+1)*lxd**3-2*toteq*ld2)
      common /ctmp1/ jduk(ld2,toteq),jdukt(toteq,ld2),scr(ctmp1pad) ! hardcoded
      integer e,eq,nstate

      real qminus(nstate,nx1,nz1,2*ndim,nelt)
      real qplus(nstate,nx1,nz1,2*ndim,nelt)
      real ujump(toteq,nx1,nz1,2*ndim,nelt)
      real nk
      integer f
      real flux(nx1*nz1)

      if (eq .eq. 1) return

      nface=2*ndim
      nxz  =nx1*nz1

      iff=1

      call rzero(flux,nxz)

      do f=1,nface
      call penalty2(qminus(1,1,1,f,e),ujump(1,1,1,f,e),nstate,flux,e,f,
     >              eq)
      call penalty2(qplus(1,1,1,f,e),ujump(1,1,1,f,e),nstate,flux,e,f,
     >              eq)
      call cmult(flux,cip(f,e),nxz)
      call copy(scr(iff),flux,nxz)
      iff=iff+nxz
      enddo ! f=1,nface

      call add_face2full_cmt(1,nx1,ny1,nz1,iface_flux(1,e),
     >                       res1(1,1,1,e,eq),scr)

      return
      end

!-----------------------------------------------------------------------

      subroutine penalty2(qface,ujump,nstate,flux,e,f,eq)
      include 'SIZE'
      include 'GEOM'
      include 'CMTDATA'

      real zenorms(lx1,lz1,6,lelt,3)
      equivalence (zenorms,unx)

      integer ctmp1pad
      parameter (ld2=lxd*lzd,ctmp1pad=(toteq+1)*lxd**3-2*toteq*ld2)
      common /ctmp1/ tmp(ld2*toteq),jdukt(toteq,ld2),notouch(ctmp1pad) ! hardcoded
      real tmp,jdukt,notouch

      integer e,f,eq ! intent(in)
      real qface(nstate,nx1,nz1)! intent(in)
      real ujump(toteq,nx1,nz1)! intent(in) 
      real flux(nx1*nz1)! intent(inout) ::
      real scr(nx1*nz1)
      real nk

      nface=2*ndim
      nxz  =nx1*nz1

      call rzero(tmp,ld2*toteq)
      call rzero(jdukt,ld2*toteq)

      il=0
      do iz=1,nz1
      do ix=1,nx1
         il=il+1
         do l=1,toteq
            jdukt(l,il)=ujump(l,ix,iz)
         enddo
      enddo
      enddo

      do j=1,ndim
         call rzero(scr,nxz)
         do k=1,ndim
! presumably, if nx1 .ne. nx1, we will have figured out how to put the appropriate
! normals in zenorms
            do i=1,nxz
               nk=zenorms(i,1,f,e,k)
               call cmult(jdukt(1,i),nk,toteq)
            enddo

! jdukt now has one face of [[Ul]]k, equation innermost. Use scr(iff)
! for visco scratch because it's huge compared to all small faces.
!        call tauij_gdu(scr,tmp,qface,qface,nstate,jdukt,jdukt,
!    >              toteq,nxz,.false.,eq,j,k,e,muf,lambdaf,tcondf,cvgf)
         enddo ! k=1,ndim

         call copy(tmp,scr,nxz)

         il=0
         do ixz=1,nxz
         il=il+1
         flux(il)=flux(il)+tmp(il)*zenorms(ixz,1,f,e,j)*area(ixz,1,f,e)
         enddo
      enddo ! j=1,3

      return
      end

! -----------------------------------------------------------------------
      subroutine compute_aux_var(e)
! this subroutine compute the surface integral contribution to auxiliary
! variable Sij=gradU and stores it in gradU. It will not be used if
! Sij=GdU unless you want gradu to be toteq times bigger
      include 'SIZE'
      include  'CMTDATA'
      include  'DG'
      include  'INPUT'
! diagnostic
      include 'GEOM'
      include 'TSTEP'
! diagnostic

      parameter (ldd=lxd*lzd*2*ldim)

      real jdu(ldd,toteq,ldim)
      integer lfq,lfc,hcsize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelcmt,
     >                   hcsize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*ldim*lfq)
      common /CMTSURFLX/ oldface(hcsize),fatface(hdsize)

      integer e ! intent(in)
      integer eq

      nfq=nx1*nz1*2*ndim*nelt
      nstate=nqq

      iqm =1
      iqp =iqm+nstate*nfq
      iuj =iqp+nstate*nfq
! get -area/bm1*[[Ul]]k refined into fine faces, jdu(:,1)
      call AuxFlux(jdu,oldface(iqm),oldface(iuj),e,nstate)
      do j=1,ndim
      do eq=1,toteq
! and add_face2full to turn totalh into Sij, i=eq
         call add_face2full_cmt(1,nx1,ny1,nz1,iface_flux(1,e),
     >                          gradu(1,1,1,eq,j),jdu(1,eq,j))
      enddo
      enddo
      return
      end
