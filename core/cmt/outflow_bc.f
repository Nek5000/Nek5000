c--------------------------------------------------------------------
      subroutine outflow(nvar,f,e,faceq,bcq,flux)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'

      integer  nvar,f,e
      real faceq(nvar,nx1,nz1,2*ldim,nelt)
      real bcq(nvar,nx1,nz1,2*ldim,nelt)
      real flux(nx1*nz1,2*ldim,nelt,*)

      call outflow_rflu(nvar,f,e,faceq,bcq,flux)

      return
      end

!--------------------------------------------------------------------

      subroutine outflow_rflu(nvar,f,e,faceq,bcq,flux1)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'CMTBCDATA'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'
      integer i,bcOpt
      integer  f,e
      real faceq(nvar,nx1*nz1,2*ldim,nelt)
      real bcq(nvar,nx1*nz1,2*ldim,nelt)
      real flux1(nx1*nz1,2*ldim,nelt,*)
      real sxn,syn,szn,rl,ul,vl,wl,pl,rr,ur,vr,wr,pr,flx(5),fs

      nxz=nx1*nz1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)     ! gives us phi- and rho-
         call userbc (ix,iy,iz,f,ieg) ! just for molarmass
         l=l+1
         sxn = unx(l,1,f,e)
         syn = uny(l,1,f,e)
         szn = unz(l,1,f,e)
         rhou= faceq(iu2,l,f,e)/phi
         rhov= faceq(iu3,l,f,e)/phi
         rhow= faceq(iu4,l,f,e)/phi
         rhoe= faceq(iu5,l,f,e)/phi
         pl  = faceq(ipr,l,f,e)
         cpgas=faceq(icpf,l,f,e)/rho
         cvgas=faceq(icvf,l,f,e)/rho
         fs = 0.0
         if(outflsub)then
            pres= pinfty
         else
            pres= pl
         endif
         call BcondOutflowPerf(1,pres,sxn,syn,szn,cpgas,molarmass,
     >                         rho,rhou,rhov,rhow,rhoe,pl,
     >                         rhob,rhoub,rhovb,rhowb,rhoeb )
         bcq(irho,l,f,e)=rhob
         bcq(iux, l,f,e)=rhoub/rhob
         bcq(iuy, l,f,e)=rhovb/rhob
         bcq(iuz, l,f,e)=rhowb/rhob
! dammit fix this. tdstate to the rescue?
         bcq(ithm,l,f,e)=(rhoeb-0.5*(rhoub**2+rhovb**2+rhowb**2)/rhob)/
     >                   cvgas
! dammit fix that
         bcq(iu1, l,f,e)=rhob*phi
         bcq(iu2, l,f,e)=rhoub*phi
         bcq(iu3, l,f,e)=rhovb*phi
         bcq(iu4, l,f,e)=rhowb*phi
         bcq(iu5, l,f,e)=rhoeb*phi
         bcq(iph, l,f,e)=phi
         bcq(ipr, l,f,e)=pres
         fs=0.0
         pl=pl*phi
         pres=pres*phi
         call CentralInviscid_FluxFunction(sxn,syn,szn,fs,
     >                                     faceq(iu1,l,f,e),pl,! U-
     >                                     bcq(iu1,l,f,e),pres,! U+,bc
     >                                     flx)
         do j=1,5
            flux1(l,f,e,j)=flx(j)
         enddo
      enddo
      enddo
      enddo

      return
      end subroutine outflow_rflu


!******************************************************************************
!
! Purpose: set outflow boundary condition for one cell.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              characteristics method of Whitfield and Janus: Three-Dimensional
!              Unsteady Euler Equations Solution Using Flux Vector Splitting.
!              AIAA Paper 84-1552, 1984. The supersonic boundary condition
!              consists of simple extrapolation.
!
! Input: bcOpt    = boundary treatment: subsonic, supersonic, or mixed
!        pout     = given static outlet pressure
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for thermally and calorically
!        perfect gas (supersonic outflow valid for all gases).
!
!******************************************************************************
!
! $Id: BcondOutflowPerf.F90,v 1.1.1.1 2014/05/05 21:47:47 tmish Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

      SUBROUTINE BcondOutflowPerf( bcOpt,pout,sxn,syn,szn,cpgas,mol,
     >                       rho,rhou,rhov,rhow,rhoe,press,
     >                       rhob,rhoub,rhovb,rhowb,rhoeb )

      IMPLICIT NONE
      integer bcopt_subsonic,bcopt_mixed
      parameter (BCOPT_SUBSONIC=1)
      parameter (BCOPT_MIXED=0)
      real MixtPerf_C_DGP,MixtPerf_Eo_DGPUVW,MixtPerf_G_CpR,MixtPerf_R_M
     >   , MixtPerf_P_DEoGVm2
      external MixtPerf_C_DGP,MixtPerf_Eo_DGPUVW,MixtPerf_G_CpR,
     >     MixtPerf_R_M,MixtPerf_P_DEoGVm2

! ... parameters
      INTEGER bcOpt

      REAL pout
      REAL rho, rhou, rhov, rhow, rhoe, press
      REAL sxn, syn, szn, cpgas, mol
      REAL rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
      REAL csound, rgas, gamma, gam1, u, v, w, mach, rrhoc, deltp,
     >            ub, vb, wb, vnd

!******************************************************************************
! gas properties; velocity components; Mach number

      rgas  = MixtPerf_R_M( mol )
      gamma = MixtPerf_G_CpR( cpgas,rgas )
      gam1  = gamma - 1.0

      u      = rhou/rho
      v      = rhov/rho
      w      = rhow/rho
      csound = MixtPerf_C_DGP( rho,gamma,press )
      mach   = SQRT(u*u+v*v+w*w)/csound

! subsonic flow ---------------------------------------------------------------

      IF (mach .lt. 1.0 .AND.
     >(bcOpt .eq. BCOPT_SUBSONIC .OR. bcOpt .eq. BCOPT_MIXED)) THEN
         rrhoc = 1.0/(rho*csound)
         deltp = press - pout
         rhob  = rho - deltp/(csound*csound)
         ub    = u + sxn*deltp*rrhoc
         vb    = v + syn*deltp*rrhoc
         wb    = w + szn*deltp*rrhoc

! - special treatment to prevent "deltp" from changing the sign
!   of velocity components. This may happen for very small u, v, w.

         vnd = ub*sxn + vb*syn + wb*szn
         IF ( vnd .lt. 0.0 ) THEN ! inflow at outflow boundary
            ub = SIGN(1.0,u)*MAX(ABS(ub),ABS(u))
            vb = SIGN(1.0,v)*MAX(ABS(vb),ABS(v))
            wb = SIGN(1.0,w)*MAX(ABS(wb),ABS(w))
         END IF ! vnd

         rhoub = rhob*ub
         rhovb = rhob*vb
         rhowb = rhob*wb
         rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pout,ub,vb,wb )

! supersonic flow -------------------------------------------------------------

      ELSE
         rhob  = rho
         rhoub = rhou
         rhovb = rhov
         rhowb = rhow
         rhoeb = rhoe
      END IF ! mach

      END SUBROUTINE BcondOutflowPerf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondOutflowPerf.F90,v $
! Revision 1.1.1.1  2014/05/05 21:47:47  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/03/26 20:21:09  haselbac
! Fix mistake in declarations
!
! Revision 1.1  2004/12/01 16:48:04  haselbac
! Initial revision after changing case
!
! Revision 1.5  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************

