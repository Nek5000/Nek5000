c-----------------------------------------------------------------------
      subroutine cmt_nek_advance
c
c     Solve the Euler equations

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'SOLN'
      include 'GEOM'
      include 'CTIMER'
      include 'CMTDATA'
      include 'CMTTIMERS'
      
      integer e,eq
      character*32 dumchars

      ftime_dum = dnekclock()
      nxyz1=lx1*ly1*lz1
      n = nxyz1*lelcmt*toteq
      nfldpart = ndim*npart

      if(istep.eq.1) call set_tstep_coef
      if(istep.eq.1) call cmt_flow_ics(ifrestart)
      if(istep.eq.1) call init_cmt_timers

      nstage = 3
      do stage=1,nstage
         if (stage.eq.1) call copy(res3(1,1,1,1,1),U(1,1,1,1,1),n)

         rhst_dum = dnekclock()
         call compute_rhs_and_dt
         rhst = rhst + dnekclock() - rhst_dum

!        if (mod(istep,res_freq).eq.0.or.istep.eq.1)then
!          dumchars='residue'
!          call dumpresidue(dumchars,stage)
!        endif
!        call exitt
c JH061114 this loop may need some work. stride difficulties

! JH111815 soon....
!        do eq=1,toteq
!           call fbinvert(res1(1,1,1,1,eq))
!        enddo

         do e=1,nelt
            do eq=1,toteq
            do i=1,nxyz1
c multiply u with bm1 as res has been multiplied by bm1 in compute_rhs
               u(i,1,1,eq,e) = bm1(i,1,1,e)*tcoef(1,stage)
     >                     *res3(i,1,1,eq,e)+bm1(i,1,1,e)*
     >                     tcoef(2,stage)*u(i,1,1,eq,e)-
     >                     tcoef(3,stage)*res1(i,1,1,e,eq)
c              u(i,1,1,eq,e) = bm1(i,1,1,e)*u(i,1,1,eq,e) - DT *
c    >                        (c1*res1(i,1,1,e,eq) + c2*res2(i,1,1,e,eq)
c    >                       + c3*res3(i,1,1,e,eq))
c-----------------------------------------------------------------------
c this completely stops working if B become nondiagonal for any reason.
! JH111815 in fact, I'd like to redo the time marching stuff above and
!          have an fbinvert call for res1
               u(i,1,1,eq,e) = u(i,1,1,eq,e)/bm1(i,1,1,e)
c that completely stops working if B become nondiagonal for any reason.
!-----------------------------------------------------------------------
            enddo
            enddo
         enddo
c        update particle rk3 stage, if there are particles
         if (arethereparticles) then
            call stokes_particles
         endif
      enddo
      call compute_primitive_vars
      ftime = ftime + dnekclock() - ftime_dum

      if (mod(istep,iostep).eq.0.or.istep.eq.1.or.istep.eq.2)then
         call out_pvar_nek
         call out_fld_nek
         call mass_balance(if3d)
      end if

      call print_cmt_timers

 101  format(4(2x,e18.9))
      return
      end

c-----------------------------------------------------------------------

      subroutine compute_rhs_and_dt
!> doxygen comments look like this
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'CMTDATA'
      include 'CTIMER'

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelcmt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*ldim*lfq)
! not sure if viscous surface fluxes can live here yet
      common /CMTSURFLX/ flux(heresize),ViscousStuff(hdsize)
      real ViscousStuff

      integer e,eq
      real wkj(lx1+lxd)
      character*32  dumchars

      if (nxd.gt.nx1) then
         call set_dealias_face
c        write(6,*)'call set dealias face'
      else
c        write(6,*)'call set alias rx'
         call set_alias_rx(istep)
      endif
c     write(6,*)'istep :', istep
!     call set_dealias_rx ! done in set_convect_cons,
! JH113015                ! now called from compute_primitive_variables

!     filter the conservative variables before start of each
!     time step
      if(IFFLTR)  call filter_cmtvar(IFCNTFILT)
!        primitive vars = rho, u, v, w, p, T, phi_g
      if (istep.eq.1) then
         call compute_primitive_vars
      else
         if(stage.gt.1) call compute_primitive_vars
      endif
!-----------------------------------------------------------------------
! JH072914 We can really only proceed with dt once we have current
!          primitive variables. Only then can we compute CFL and/or dt.
!          I know this isn't an ideal place for it, but it avoids some
!          repeated work. If we ever go to RK, an if(isubstage==1) will
!          be needed here.
!-----------------------------------------------------------------------
      if(stage.eq.1) call setdtcmt

!     !Total_eqs = 5 (we will set this up so that it can be a user 
!     !defined value. 5 will be its default value)
!     !eq = 1 -------- Mass balance
!     !eq = 2 -------- x  momentum 
!     !eq = 3 -------- y  momentum 
!     !eq = 4 -------- z  momentum 
!     !eq = 5 -------- Energy Equation 

!-----------------------------------------------------------------------
! JH060314 Compute face fluxes now that we have the primitive
!          variables. As face selection gets smarter (e.g.: interior
!          faces only), surface_fluxes_inviscid's argument list may grow
!-----------------------------------------------------------------------
      call fluxes_full_field
      ntot = lx1*ly1*lz1*lelcmt*toteq
      call rzero(res1,ntot)

      nstate=nqq
      nfq=nx1*nz1*2*ndim*nelt
      iqm =1
      iqp =iqm+nstate*nfq
      iflx=iqp+nstate*nfq
      do eq=1,toteq
         ieq=(eq-1)*ndg_face+iflx
         call surface_integral_full(res1(1,1,1,1,eq),flux(ieq))
      enddo
      iuj=iflx ! overwritten with [[U]]

!     if (ifvisc) call ujump ! need to write something that goes from
!     U+- in fatface to ujump
      if (ifvisc) call compute_transport_props

      do e=1,nelt
! Get user defined forcing from userf defined in usr file
         call cmtusrf(e)
         if (ifvisc)then
             call compute_gradients(e)
! compute_aux_var will likely not be called if Sij=GdU
             call compute_aux_var(e)
!------------------------------
! NB!!!!! gradu=Sij, NOT dUl/dxk!!!!
!------------------------------
         endif
         do eq=1,toteq
            call assemble_h(e,eq)
! compute the volume integral term and add to res1(:,e,eq)
            if (nxd.gt.nx1) then
               call flux_div_integral_dealiased(e,eq)
            else
               call flux_div_integral_aliased(e,eq)
            endif
!------------------------------
! JH050615 BR1 ONLY for now
!           if (.not.ifbr1)
!    >      call penalty(flux(iqm),flux(iqp),flux(iuj),e,eq,nstate)
!------------------------------
! Compute the forcing term in each of the 5 eqs
            call compute_forcing(e,eq)
         enddo
      enddo

! get the other half of Hij^{d*}
      if (ifvisc)then
         call viscousf
         do eq=2,toteq ! no viscous flux of mass
            ieq=(eq-1)*ndg_face+1
!Finally add viscous surface flux functions of derivatives to res1.
            call surface_integral_full(res1(1,1,1,1,eq),flux(ieq))
         enddo
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine set_tstep_coef
      include 'SIZE'
      include 'TSTEP'

      COMMON /TIMESTEPCOEf/ tcoef(3,3)

      COMMON /TSTEPSTAGE/  stage,nstage
      integer              stage,nstage

      nstage = 3
      do stage=1,nstage
         if(stage.eq.1)then
            tcoef(1,stage) = 0.0
            tcoef(2,stage) = 1.0 
            tcoef(3,stage) = dt 
         endif
         if(stage.eq.2)then
            tcoef(1,stage) = 3.0/4.0
            tcoef(2,stage) = 1.0/4.0 
            tcoef(3,stage) = dt/4.0 
         endif
         if(stage.eq.3)then
            tcoef(1,stage) = 1.0/3.0
            tcoef(2,stage) = 2.0/3.0 
            tcoef(3,stage) = dt*2.0/3.0 
         endif
      enddo
      return
      end
!-----------------------------------------------------------------------

      subroutine cmt_flow_ics(ifrestart)
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'

      logical ifrestart
      integer e
      nxyz1 = nx1*ny1*nz1
      n     = nxyz1*lelcmt*toteq
      if (ifrestart)then
         do e=1,nelt
            call copy(U(1,1,1,2,e),vx(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,3,e),vy(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,4,e),vz(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,5,e),t(1,1,1,e,1),nxyz1) 
            call copy(U(1,1,1,1,e),pr(1,1,1,e),nxyz1) 
         enddo
      endif
      call rzero(res1,n)
      call rzero(res2,n)
      return
      end
!-----------------------------------------------------------------------

      subroutine print_cmt_timers
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'

c we need our own IO features. Until then we use the default nek routines
      if ((mod(istep,flio_freq).eq.0.and.istep.gt.0)
     $                               .or.istep.eq.nstep)then
         dmtime1 = ftime/istep
         dtime_ = glsum(dmtime1,1)
         if(nio.eq.0) write(6,*) 'fluid rhs compute time(Avg)  '
     $               ,dtime_/np
      endif
      return 
      end
!-----------------------------------------------------------------------

      subroutine init_cmt_timers
      include 'CMTTIMERS'

      rhst    = 0.00
      ftime   = 0.00

      return
      end
!-----------------------------------------------------------------------
