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

      call InviscidFlux(fatface(iqm),fatface(iqp),fatface(iflx)
     >                 ,nstate,toteq)

!     call face_flux_commo(fatface(iflx),fatface(iflx),ndg_face,toteq,
!    >                     flux_hndl) ! for non-symmetric gs_op someday

      return
      end

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
      end

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
         qminus(i,ivar)=yourface(i,1,1,1)
      enddo

      return
      end

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
      call gs_op_fields(handle,yours,nf,nstate,1,1,0)
      call sub2 (yours,mine,ntot)
      return
      end

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
      end

!-------------------------------------------------------------------------------

      subroutine InviscidFlux(qminus,qplus,flux,nstate,nflux)
!-------------------------------------------------------------------------------
! JH091514 A fading copy of RFLU_ModAUSM.F90 from RocFlu
!-------------------------------------------------------------------------------

!#ifdef SPEC
!      USE ModSpecies, ONLY: t_spec_type
!#endif
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real qminus(nx1*nz1,2*ndim,nelt,nstate),
     >     qplus(nx1*nz1,2*ndim,nelt,nstate),
     >     flux(nx1*nz1,2*ndim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface,ifield
      parameter (lfd=lxd*lzd)
! JH111815 legacy rocflu names.
!
! nx,ny,nz : outward facing unit normal components
! fs       : face speed. zero until we have moving grid
! jaco_c   : fdim-D GLL grid Jacobian
! nm       : jaco_c, fine grid
!
! State on the interior (-, "left") side of the face
! rl       : density
! ul,vl,wl : velocity
! tl       : temperature
! al       : sound speed
! pl       : pressure, then phi
! cpl      : rho*cp
! State on the exterior (+, "right") side of the face
! rr       : density
! ur,vr,wr : velocity
! tr       : temperature
! ar       : sound speed
! pr       : pressure
! cpr      : rho*cp

      COMMON /SCRNS/ nx(lfd), ny(lfd), nz(lfd), rl(lfd), ul(lfd),
     >               vl(lfd), wl(lfd), pl(lfd), tl(lfd), al(lfd),
     >               cpl(lfd),rr(lfd), ur(lfd), vr(lfd), wr(lfd),
     >               pr(lfd),tr(lfd), ar(lfd),cpr(lfd), fs(lfd),
     >               jaco_f(lfd),flx(lfd,toteq),jaco_c(lx1*lz1)
      real nx, ny, nz, rl, ul, vl, wl, pl, tl, al, cpl, rr, ur, vr, wr,
     >                pr,tr, ar,cpr, fs,jaco_f,flx,jaco_c

!     REAL vf(3)
      real nTol
      character*132 deathmessage
      character*3 cb

      nTol = 1.0E-14

      fdim=ndim-1
      nface = 2*ndim
      nxz   = nx1*nz1
      nxzd  = nxd*nzd
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

            do j=1,nstate
               do i=1,nxz
                  if (abs(qplus(i,f,e,j)) .gt. ntol) then
                  write(6,*) nid,j,i,qplus(i,f,e,j),qminus(i,f,e,j),cb,
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
            if (cb.eq.'v  ' .or. cb .eq. 'V  ') then
              call inflow(nstate,f,e,qminus,qplus,flux)
            elseif (cb.eq.'O  ') then
              call outflow(nstate,f,e,qminus,qplus,flux)
            elseif (cb .eq. 'W  ' .or. cb .eq.'I  '.or.cb .eq.'SYM')then
              call wallbc(nstate,f,e,qminus,qplus,flux)
            endif 

         else ! cbc(f,e,ifield) == 'E  ' or 'P  ' below; interior face

! JH111715 now with dealiased surface integrals. I am too lazy to write
!          something better
            call map_faced(nx,unx(1,1,f,e),nx1,nxd,fdim,0)
            call map_faced(ny,uny(1,1,f,e),nx1,nxd,fdim,0)
            call map_faced(nz,unz(1,1,f,e),nx1,nxd,fdim,0)

            call map_faced(rl,qminus(1,f,e,irho),nx1,nxd,fdim,0)
            call map_faced(ul,qminus(1,f,e,iux),nx1,nxd,fdim,0)
            call map_faced(vl,qminus(1,f,e,iuy),nx1,nxd,fdim,0)
            call map_faced(wl,qminus(1,f,e,iuz),nx1,nxd,fdim,0)
            call map_faced(pl,qminus(1,f,e,ipr),nx1,nxd,fdim,0)
            call map_faced(tl,qminus(1,f,e,ithm),nx1,nxd,fdim,0)
            call map_faced(al,qminus(1,f,e,isnd),nx1,nxd,fdim,0)
            call map_faced(cpl,qminus(1,f,e,icpf),nx1,nxd,fdim,0)

            call map_faced(rr,qplus(1,f,e,irho),nx1,nxd,fdim,0)
            call map_faced(ur,qplus(1,f,e,iux),nx1,nxd,fdim,0)
            call map_faced(vr,qplus(1,f,e,iuy),nx1,nxd,fdim,0)
            call map_faced(wr,qplus(1,f,e,iuz),nx1,nxd,fdim,0)
            call map_faced(pr,qplus(1,f,e,ipr),nx1,nxd,fdim,0)
            call map_faced(tr,qplus(1,f,e,ithm),nx1,nxd,fdim,0)
            call map_faced(ar,qplus(1,f,e,isnd),nx1,nxd,fdim,0)
            call map_faced(cpr,qplus(1,f,e,icpf),nx1,nxd,fdim,0)

            call invcol3(jaco_c,area(1,1,f,e),wghtc,nxz)
            call map_faced(jaco_f,jaco_c,nx1,nxd,fdim,0) 
            call col2(jaco_f,wghtf,nxzd)
            call rzero(fs,nxzd) ! moving grid stuff later

            call AUSM_FluxFunction(nxzd,nx,ny,nz,jaco_f,fs,rl,ul,vl,wl,
     >                        pl,al,tl,rr,ur,vr,wr,pr,ar,tr,flx,cpl,cpr)

            call map_faced(pl,qminus(1,f,e,iph),nx1,nxd,fdim,0)
            do j=1,toteq
               call col2(flx(1,j),pl,nxzd)
               call map_faced(flux(1,f,e,j),flx(1,j),nx1,nxd,fdim,1)
            enddo

         endif ! cbc(f,e,ifield)
      enddo
      enddo

      end

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

! weak form until we get the time loop rewritten
!     onem=-1.0
!     ntot=nx1*nz1*2*ndim*nelt
!     call cmult(flux,onem,ntot)
! weak form until we get the time loop rewritten
      call add_face2full_cmt(nelt,nx1,ny1,nz1,iface_flux,vol,flux)

      return
      end

!-------------------------------------------------------------------------------

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
      include  'GEOM'
      integer e,nstate
      real hsijstar(nxd,nzd,2*ndim,toteq,ldim) ! starts life as scratch
      real qminus(nstate,nx1,nz1,2*ndim,nelt),
     >        ujump(toteq,nx1,nz1,2*ndim,nelt),
     >        scrth1(lx1,lz1,2*ldim)
! JH050615 This subroutine's arguments must change if Sij=GdU
!     real qminus(nstate,nx1,nz1,2*ndim,nelt),
!    >        ujump(toteq,nx1,nz1,2*ndim,nelt),scrf1(nx1,nz1,2*ndim),
!    >        gdu(toteq,nx1,nz1)

      real zenorms(lx1,lz1,6,lelt,3)
      equivalence (zenorms,unx)

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
      end

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
      end

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
      end

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
      end
