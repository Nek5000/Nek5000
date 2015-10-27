! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities ! NOT USED IN CMT-NEK YET
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

      SUBROUTINE AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,
     >                                Hl,al,rr,ur,vr,wr,pr,Hr,ar,flx,vf)

!     IMPLICIT NONE ! HAHAHHAHHAHA

! ==============================================================================
! Arguments
! ==============================================================================

      REAL al,ar,fs,Hl,Hr,nm,nx,ny,nz,pl,pr,rl,rr,ul,ur,! INTENT(IN) ::
     >                    vl,vr,wl,wr
      REAL flx(5),vf(3) ! INTENT(OUT) ::

! ==============================================================================
! Locals
! ==============================================================================

      REAL af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr,
     >        wtl,wtr

! ******************************************************************************
! Start, compute face state
! ******************************************************************************

      ql = ul*nx + vl*ny + wl*nz - fs
      qr = ur*nx + vr*ny + wr*nz - fs

      af = 0.5*(al+ar) ! NOTE not using original formulation, see note

      ml  = ql/af
      mla = ABS(ml)

      mr  = qr/af
      mra = ABS(mr)    

      IF ( mla <= 1.0 ) THEN 
         mlp = 0.25*(ml+1.0)*(ml+1.0) + 0.125*(ml*ml-1.0)*(ml*ml-1.0)
         wtl = 0.25*(ml+1.0)*(ml+1.0)*(2.0-ml) +
     >         0.1875*ml*(ml*ml-1.0)*(ml*ml-1.0)
      ELSE
         mlp = 0.5*(ml+mla)
         wtl = 0.5*(1.0+ml/mla)
      END IF ! mla

      IF ( mra <= 1.0 ) THEN 
         mrm = -0.25*(mr-1.0)*(mr-1.0)-0.125*(mr*mr-1.0)*(mr*mr-1.0)
         wtr = 0.25*(mr-1.0)*(mr-1.0)*(2.0+mr) -
     >         0.1875*mr*(mr*mr-1.0)*(mr*mr-1.0)
      ELSE
         mrm = 0.5*(mr-mra)
         wtr = 0.5*(1.0-mr/mra)
      END IF ! mla

      mf  = mlp + mrm
      mfa = ABS(mf)
      mfp = 0.5*(mf+mfa)
      mfm = 0.5*(mf-mfa)

      pf = wtl*pl + wtr*pr 

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

      vf(1) = mfp*ul + mfm*ur
      vf(2) = mfp*vl + mfm*vr
      vf(3) = mfp*wl + mfm*wr    

      flx(1) = (af*(mfp*rl    + mfm*rr   )        )*nm
      flx(2) = (af*(mfp*rl*ul + mfm*rr*ur) + pf*nx)*nm
      flx(3) = (af*(mfp*rl*vl + mfm*rr*vr) + pf*ny)*nm
      flx(4) = (af*(mfp*rl*wl + mfm*rr*wr) + pf*nz)*nm
      flx(5) = (af*(mfp*rl*Hl + mfm*rr*Hr) + pf*fs)*nm

! ******************************************************************************
! End
! ******************************************************************************
      return
      END SUBROUTINE AUSM_FluxFunction

!-----------------------------------------------------------------------

      SUBROUTINE CentralInviscid_FluxFunction(nx,ny,nz,fs,ul,pl,
     >                                     ur,pr,flx)
! JH081915 More general, more obvious
      real nx,ny,nz,fs,ul(5),pl,ur(5),pr ! intent(in)
      real flx(5)! intent(out),dimension(5) ::

      rl =ul(1)
      rul=ul(2)
      rvl=ul(3)
      rwl=ul(4)
      rel=ul(5)

      rr =ur(1)
      rur=ur(2)
      rvr=ur(3)
      rwr=ur(4)
      rer=ur(5)

      ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
      qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

      flx(1) = 0.5*(ql* rl                + qr* rr               )
      flx(2) = 0.5*(ql* rul       + pl*nx + qr* rur       + pr*nx)
      flx(3) = 0.5*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)
      flx(4) = 0.5*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)
      flx(5) = 0.5*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)

      return
      end subroutine CentralInviscid_FluxFunction
