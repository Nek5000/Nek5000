      subroutine fluxes_full_field
!-----------------------------------------------------------------------
! JH060314 First, compute face fluxes now that we have the primitive variables
! JH091514 renamed from "surface_fluxes_inviscid" since it handles all fluxes
!          that we compute from variables stored for the whole field (as
!          opposed to one element at a time).
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelcmt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*ldim*lfq)
! JH070214 OK getting different answers whether or not the variables are
!          declared locally or in common blocks. switching to a different
!          method of memory management that is more transparent to me.
      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      integer eq
      character*32 cname
      nfq=nx1*nz1*2*ndim*nelt
      nstate = nqq
      ntot1 = nfq*nstate
! where different things live
      iqm =1
      iqp =iqm+nstate*nfq
      iflx=iqp+nstate*nfq

      call fillq(irho,vtrans,fatface(iqm),fatface(iqp))
      call fillq(iux, vx,    fatface(iqm),fatface(iqp))
      call fillq(iuy, vy,    fatface(iqm),fatface(iqp))
      call fillq(iuz, vz,    fatface(iqm),fatface(iqp))
      call fillq(ipr, pr,    fatface(iqm),fatface(iqp))
      call fillq(ithm,t,     fatface(iqm),fatface(iqp))
      call fillq(isnd,csound,fatface(iqm),fatface(iqp))
      call fillq(iph, phig,  fatface(iqm),fatface(iqp))
      call fillq(icvf,vtrans(1,1,1,1,icv),fatface(iqm),fatface(iqp))
      call fillq(icpf,vtrans(1,1,1,1,icp),fatface(iqm),fatface(iqp))
      if (ifvisc) then
         call fillq(imuf, vdiff(1,1,1,1,imu), fatface(iqm),fatface(iqp))
         call fillq(ikndf,vdiff(1,1,1,1,iknd),fatface(iqm),fatface(iqp))
         call fillq(ilamf,vdiff(1,1,1,1,ilam),fatface(iqm),fatface(iqp))
      endif

      i_cvars=(iu1-1)*nfq+1
      do eq=1,toteq
         call faceu(eq,fatface(i_cvars))
         i_cvars=i_cvars+nfq
      enddo

      call face_state_commo(fatface(iqm),fatface(iqp),nfq,nstate
     >                     ,dg_hndl)

! JH091614 Since we don't have a correctly vectorized gs_op_fields here,
!          I have decided to transpose the q's to get primitive vars
!          ordered innermost.
      call transpose(fatface(iflx),nstate,fatface(iqm),nfq)
      call copy(fatface(iqm),fatface(iflx),ntot1)
      call transpose(fatface(iflx),nstate,fatface(iqp),nfq)
      call copy(fatface(iqp),fatface(iflx),ntot1)

      call InviscidFlux(fatface(iqm),fatface(iqp),fatface(iflx)
     >                 ,nstate,toteq)

!     call face_flux_commo(fatface(iflx),fatface(iflx),ndg_face,toteq,
!    >                     flux_hndl) ! for non-symmetric gs_op someday

      return
      end subroutine fluxes_full_field

!-----------------------------------------------------------------------

      subroutine faceu(ivar,yourface)
! get faces of conserved variables stored contiguously
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'
      integer e
      real yourface(nx1,nz1,2*ldim,nelt)

      do e=1,nelt
         call full2face_cmt(1,nx1,ny1,nz1,iface_flux(1,e),
     >                      yourface(1,1,1,e),u(1,1,1,ivar,e))
      enddo

      return
      end subroutine faceu

!-----------------------------------------------------------------------

      subroutine fillq(ivar,field,qminus,yourface)
      include 'SIZE'
      include 'DG'

      integer ivar! intent(in)
      real field(nx1,ny1,nz1,nelt)! intent(in)
!     real, intent(out)qminus(7,nx1*nz1*2*ldim*nelt) ! gs_op no worky
      real qminus(nx1*nz1*2*ndim*nelt,*)! intent(out)
      real yourface(nx1,nz1,2*ndim,*)
      integer e,f

      nxz  =nx1*nz1
      nface=2*ndim

      call full2face_cmt(nelt,nx1,ny1,nz1,iface_flux,yourface,field)

      do i=1,ndg_face
!        qminus(ivar,i)=yourface(i,1,1,1) ! gs_op_fields no worky yet.
                                          ! tranpose later
         qminus(i,ivar)=yourface(i,1,1,1)
      enddo

      return
      end subroutine fillq

!-----------------------------------------------------------------------

      subroutine face_state_commo(mine,yours,nf,nstate,handle)

! JH060414 if we ever want to be more intelligent about who gets what,
!          who gives what and who does what, this is the place where all
!          that is done. At the very least, gs_op may need the transpose
!          flag set to 1. Who knows. Everybody duplicates everything for
!          now.
! JH070714 figure out gs_op_fields, many, vec, whatever (and the
!          corresponding setup) to get this done for the transposed
!          ordering of state variables. I want variable innermost, not
!          grid point.

      integer handle,nf,nstate ! intent(in)
      real yours(*),mine(*)

      ntot=nf*nstate
      call copy(yours,mine,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
!     write(6,*) 'face_state_commo ',nstate
      call gs_op_fields(handle,yours,nf,nstate,1,1,0)
      call sub2 (yours,mine,ntot)
      return
      end subroutine face_state_commo

!-----------------------------------------------------------------------

      subroutine face_flux_commo(flux1,flux2,nf,neq,handle)
! JH060514 asymmetric transposed gs_op, gs_unique magic may be needed if
!          we ever decide to avoid redundancy. For now, this routine
!          doesn't need to do anything.
      integer ntot,handle
      real flux1(*),flux2(*)
! JH061814 It doesn't need to do anything, but a sanity check would be
!          wise.
      return
      end subroutine face_flux_commo

!-------------------------------------------------------------------------------

      subroutine InviscidFlux(qminus,qplus,flux,nstate,nflux)
!-------------------------------------------------------------------------------
! JH091514 A fading copy of RFLU_ModAUSM.F90 from RocFlu
!-------------------------------------------------------------------------------

!#ifdef SPEC
!      USE ModSpecies, ONLY: t_spec_type
!#endif
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'CMTDATA'
      include 'DG'
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************
      real MixtPerf_Ho_CpTUVW
      external MixtPerf_Ho_CpTUVW

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real qminus(nstate,nx1*nz1,2*ndim,nelt),
     >     qplus(nstate,nx1*nz1,2*ndim,nelt),
     >     flux(nx1*nz1,2*ndim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,i,k,nxz,nface,ivoid,ifield
      REAL al,ar,cpl,cpr,dx,dy,dz,fs,fsu,gcl,gcr,gl,gr,Hl,Hr,irl,irr,
     >      nm,nTol,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,rgl,rgr
      REAL flx(5),vf(3)
      character*132 deathmessage
      character*3 cb

      nTol = 1.0E-14

      nface = 2*ndim
      nxz   = nx1*nz1
      ifield= 1

!     if (outflsub)then
!        call maxMachnumber
!     endif
      do e=1,nelt
      do f=1,nface

         cb=cbc(f,e,ifield)
         if (cb.ne.'E  '.and.cb.ne.'P  ') then ! cbc bndy

!-----------------------------------------------------------------------
! compute flux for weakly-enforced boundary condition
!-----------------------------------------------------------------------

            do i=1,nxz
               do j=1,nstate
                  if (abs(qplus(j,i,f,e)) .gt. ntol) then
                  write(6,*) nid,j,i,qplus(j,i,f,e),qminus(j,i,f,e),cb,
     > nstate
                  write(deathmessage,*)  'GS hit a bndy,f,e=',f,e,'$'
! Make sure you are not abusing this error handler
                  call exitti(deathmessage,f)
                  endif
               enddo
            enddo
! JH031315 flux added to argument list. BC routines preserve qminus for
!          obvious reasons and fill qplus with good stuff for everybody:
!          imposed states for Dirichlet conditions, and important things
!          for viscous numerical fluxes.
! JH060215 added SYM bc. Just use it as a slip wall hopefully.
            if (cb.eq.'v  ') then
              call inflow(nstate,f,e,qminus,qplus,flux)
            elseif (cb.eq.'O  ') then
              call outflow(nstate,f,e,qminus,qplus,flux)
            elseif (cb .eq. 'W  ' .or. cb .eq.'I  '.or.cb .eq.'SYM')then
              call wallbc(nstate,f,e,qminus,qplus,flux)
            endif 

         else ! cbc(f,e,ifield) == 'E  ' or 'P  ' below; interior face

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

            do i=1,nxz
               do j=1,5
                  flux(i,f,e,j)=0.0
               enddo

! ==============================================================================
!   Get face geometry and, when it exists, grid speed
! ==============================================================================

               nx = unx(i,1,f,e)
               ny = uny(i,1,f,e)
               nz = unz(i,1,f,e)
               nm = 1.0 ! multiply by area(i,1,f,e) later

!              fs = pGrid%gs(indGs*ifg) ! face speed
               fs = 0.0 ! moving grid stuff later

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

               rl = qminus(irho,i,f,e)
               ul = qminus(iux,i,f,e)
               vl = qminus(iuy,i,f,e)
               wl = qminus(iuz,i,f,e)
               pl = qminus(ipr,i,f,e)
               tl = qminus(ithm,i,f,e)
               al = qminus(isnd,i,f,e)
               cpl= qminus(icpf,i,f,e)/rl
               Hl = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

               rr = qplus(irho,i,f,e)
               ur = qplus(iux,i,f,e)
               vr = qplus(iuy,i,f,e)
               wr = qplus(iuz,i,f,e)
               pr = qplus(ipr,i,f,e)
               tr = qplus(ithm,i,f,e)
               ar = qplus(isnd,i,f,e)
               cpr= qplus(icpf,i,f,e)/rr
               Hr = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)

! ==============================================================================
!   Compute fluxes
! ==============================================================================

               call AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,
     >                                   al,rr,ur,vr,wr,pr,Hr,ar,flx,vf)
               do j=1,5
                  flux(i,f,e,j)=flx(j) ! put phi here soon
               enddo

            enddo ! i=1,nxz

         endif ! cbc(f,e,ifield)
      enddo
      enddo

      END subroutine InviscidFlux

!-----------------------------------------------------------------------

      subroutine surface_integral_full(vol,flux)
! Integrate surface fluxes for an entire field. Add contribution of flux
! to volume according to add_face2full_cmt
      include 'SIZE'
      include 'GEOM'
      include 'DG'
      include 'CMTDATA'
      real vol(nx1*ny1*nz1*nelt),flux(*)
      integer e,f

      nxz=nx1*nz1
      nface=2*ldim
      k=0

! JH030915 As qminus and qplus grow and get more flexible, throttling by phig
!          should be done via qminus and NOT the way it is done here. We don't
!          need icpvars anymore, and, most importantly ViscousFlux ALREADY HAS PHI!!!!!!
      l = 0
      do e=1,nelt
      do f=1,nface
         call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f) 
         k = 0
         do iz=k0,k1
         do iy=j0,j1
         do ix=i0,i1
            l = l + 1
            k = k + 1
            flux(l)=flux(l)*area(k,1,f,e)*phig(ix,iy,iz,e)
         enddo
         enddo
         enddo
      enddo
      enddo
      call add_face2full_cmt(nelt,nx1,ny1,nz1,iface_flux,vol,flux)

      return
      end subroutine surface_integral_full

!-------------------------------------------------------------------------------

      subroutine surface_integral_elm(e,eq)
! JH062314 surface integrals one element at a time. Only needed because we are
!          doing all this in strong form
! JH070214 CMTSURFLX no longer used here. it's just one element, we can afford
!          to store stuff that lives here only
! JH091015 this routine had better vanish when we get weak form working
      include 'SIZE'
      include 'DG'
      include 'CMTDATA'
      include 'INPUT'

      common /GSURF/ ZENORMS(LX1,LZ1,6,lelcmt,3)
     >             ,SUPAHUGET(LX1*LZ1*6*lelcmt*6)
     >             ,AREA  (LX1,LZ1,6,lelcmt)
     >             ,DLAM

      integer lfc1
      parameter (lfc1=lx1*lz1*2*ldim)
      real flux(lfc1),yrface(lfc1)
      integer e,eq
      integer f

      nface = 2*ndim
      nxz   = nx1*nz1

      call rzero(flux,lfc1)
      do j=1,ndim
         call full2face_cmt(1,nx1,ny1,nz1,iface_flux(1,e),yrface,
     >                      totalh(1,j))
! -n-.H
         k=0
         do f=1,nface
            do i=1,nxz
            k=k+1
            flux(k) = flux(k)-area(i,1,f,e)*zenorms(i,1,f,e,j)*yrface(k)
            enddo
         enddo

      enddo

      call add_face2full_cmt(1,nx1,ny1,nz1,iface_flux(1,e),
     >                       res1(1,1,1,e,eq),flux)

      return
      end subroutine surface_integral_elm

!-----------------------------------------------------------------------

! JH050615 This subroutine's arguments must change if Sij=GdU
!     subroutine AuxFlux(hsijstar,qminus,ujump,scrf1,gdu,e,eq,j,
!    >                   nstate)
      subroutine AuxFlux(hsijstar,qminus,ujump,e,nstate)
! MS050615 Computes jump penalization on element faces for auxiliary
!          variable Sij=gradU
      include  'SIZE'
      include  'MASS'
      include  'DG'
      include  'CMTDATA'
      include  'INPUT'
      integer e,nstate
      real hsijstar(nxd,nzd,2*ndim,toteq,ldim) ! starts life as scratch
      real qminus(nstate,nx1,nz1,2*ndim,nelt),
     >        ujump(toteq,nx1,nz1,2*ndim,nelt),
     >        scrth1(lx1,lz1,2*ldim)
! JH050615 This subroutine's arguments must change if Sij=GdU
!     real qminus(nstate,nx1,nz1,2*ndim,nelt),
!    >        ujump(toteq,nx1,nz1,2*ndim,nelt),scrf1(nx1,nz1,2*ndim),
!    >        gdu(toteq,nx1,nz1)
      common /GSURF/ ZENORMS(LX1,LZ1,6,lelcmt,3)
     >             ,SUPAHUGET(LX1*LZ1*6*lelcmt*6)
     >             ,AREA  (LX1,LZ1,6,lelcmt)
     >             ,DLAM
      integer f

      nface=2*ndim
      nxz  =nx1*nz1
      ntot =nxz*nface

      do k=1,ndim
         do l=1,toteq
         do f=1,nface
            do i=1,nxz
               hsijstar(i,1,f,l,k)=-0.5*ujump(l,i,1,f,e)*area(i,1,f,e)
     >                                               *zenorms(i,1,f,e,k)
            enddo
         enddo
         enddo
      enddo

      call full2face_cmt(1,nx1,ny1,nz1,iface_flux(1,e),scrth1(1,1,1),
     >                   bm1(1,1,1,e))
      do k=1,ndim
         do l=1,toteq
            call invcol2(hsijstar(1,1,1,l,k),scrth1(1,1,1),ntot)
         enddo 
      enddo
      return
      end subroutine AuxFlux

!-----------------------------------------------------------------------

      subroutine store_gdu_hstarsij(gdudxk,hstarsij,e,eq,j)
      include 'SIZE'
      include 'DG'
      include 'CMTDATA'
      include 'INPUT'
      integer e,eq,j
      real gdudxk(nx1,nz1,2*ndim,nelt,toteq,ndim)
      real hstarsij(nx1,ny1,nz1)
      integer f

      nface = 2*ndim
      nxz   = nx1*nz1

      call full2face_cmt(1,nxd,nyd,nzd,iface_flux(1,e),
     >                   gdudxk(1,1,1,e,eq,j),hstarsij)
      
      call cmult(gdudxk(1,1,1,e,eq,j),-1.0,nxz*nface)

      return
      end subroutine store_gdu_hstarsij

!-----------------------------------------------------------------------

      subroutine viscousf
! gets central-flux contribution to interior penalty numerical flux
! Hij^{d*}
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'

      integer lfq,lfc,hcsize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelcmt,
     >                   hcsize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*ldim*lfq)
      common /CMTSURFLX/ flux(hcsize),gdudxk(hdsize)
      real flux
      real const
      integer e,f

      nface=2*ndim
      nfc =nx1*nz1*nface*nelt
      nnfc=toteq*ndim
      ntot=nfc*nnfc
      const = 0.5
      nxz = nx1*nz1
      call cmult(gdudxk,const,ntot)
!-----------------------------------------------------------------------
! supa huge gs_op to get {{gdu}}
! operation flag is second-to-last arg, an integer
!                                                   1 ==> +
      call gs_op_fields(dg_hndl,gdudxk,nfc,nnfc,1,1,0)
!-----------------------------------------------------------------------

! now dot with n-
      call special_dot_n(flux,gdudxk)
      call bcflux(flux) ! needs work

      return
      end subroutine viscousf

!-----------------------------------------------------------------------

      subroutine special_dot_n(fdotn,face)
! no, I don't feel like making this routine more general
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      real fdotn(nx1*nz1*2*ndim*nelt,toteq)
      real face (nx1*nz1*2*ndim*nelt,toteq,ndim) ! must be 3?
      integer e,eq

      nxz=nx1*nz1
      nface=2*ndim
      ntot=nxz*nface

      do eq=1,toteq
         l=1
         do e=1,nelt
            call col3   (fdotn(l,eq),face(l,eq,1),unx(1,1,1,e),ntot)
            l=l+ntot
         enddo
      enddo

      do eq=1,toteq
         l=1
         do e=1,nelt
            call addcol3(fdotn(l,eq),face(l,eq,2),uny(1,1,1,e),ntot)
            l=l+ntot
         enddo
      enddo

      if (if3d) then
         do eq=1,toteq
            l=1
            do e=1,nelt
               call addcol3(fdotn(l,eq),face(l,eq,3),unz(1,1,1,e),ntot)
               l=l+ntot
            enddo
         enddo
      endif

      return
      end subroutine special_dot_n
