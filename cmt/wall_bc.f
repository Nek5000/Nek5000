      subroutine wallbc(nstate,f,e,faceq,bcq,flux)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      include 'CMTBCDATA'

      integer nstate,f,e
      real    faceq(nstate,nx1*nz1,2*ndim,nelt)
      real    bcq(nstate,nx1*nz1,2*ndim,nelt) 
      real    flux(nx1*nz1,2*ndim,nelt,*)

      call slipwall_rflu(nstate,f,e,faceq,bcq,flux) ! calls RindState stuff

      return
      end subroutine wallbc

!--------------------------------------------------------------------

      subroutine slipwall_rflu(nvar,f,e,faceq,bcq,fluxw)
      include 'SIZE'
      include 'CMTBCDATA'
      include 'GEOM'
      include 'NEKUSE'
      include 'INPUT'
      include 'PARALLEL'
      integer  f,e
! JH091614 faceq now has intent(inout)...
! JH031315 not anymore. nobody changes qminus here. that's dumb
      real faceq(nvar,nx1*nz1,2*ndim,nelt)
      real bcq(nvar,nx1*nz1,2*ndim,nelt)
      real fluxw(nx1*nz1,2*ndim,nelt,*)
      integer i, nxz
      real nx,ny,nz,rl,ul,vl,wl,pl,fs,dumminus(5),flx(5)

      nxz=nx1*nz1
      ieg=lglel(e)

! I know this says slipwall, but to the inviscid terms all walls are
! slip. or something.
      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg)
         l=l+1
         nx = unx(l,1,f,e)
         ny = uny(l,1,f,e)
         nz = unz(l,1,f,e)
         rl = faceq(irho,l,f,e)
         ul = faceq(iux,l,f,e)
         vl = faceq(iuy,l,f,e)
         wl = faceq(iuz,l,f,e)
         pl = faceq(ipr,l,f,e)
         fs = 0.0 ! no moving grid for awhile, and it will not look anything
                  ! like RocFlu
         call RFLU_SetRindStateSlipWallPerf(cp,molarmass,nx,ny,nz,
     >                                      rl,ul,vl,wl,fs,pl)
         bcq(irho,l,f,e)=rl
         bcq(iux,l,f,e)=ul
         bcq(iuy,l,f,e)=vl
         bcq(iuz,l,f,e)=wl
         bcq(ipr,l,f,e)=pl ! from RFLU_SetRindStateSlipWallPerf
         bcq(iph,l,f,e)=phi
         if (ifvisc) then
            bcq(iux,l,f,e)=ux
            bcq(iuy,l,f,e)=uy
            bcq(iuz,l,f,e)=uz
            if (cbc(f,e,2) .eq. 'W  ') bcq(ithm,l,f,e)=temp
         endif
! Inviscid flux at walls is due to pressure only. should probably just
! hardcode that instead of calling CentralInviscid so trivially
         dumminus(1)=rl
         dumminus(2)=0.0
         dumminus(3)=0.0
         dumminus(4)=0.0
         dumminus(5)=0.0
         call CentralInviscid_FluxFunction(nx,ny,nz,fs,dumminus,pl
     >                                    ,dumminus,pl,flx)
         do j=1,5
            fluxw(l,f,e,j)=flx(j)
         enddo
         if (cbc(f,e,2) .eq. 'I  ') fluxw(l,f,e,5) = flux
      enddo
      enddo
      enddo

      return
      end subroutine slipwall_rflu

!-----------------------------------------------------------------------
! ******************************************************************************
!
! Purpose: Set rind state for slip-wall boundaries and perfect gas.
!
! Description: Torn bleeding from RocFlu. I think "rind" means the same thing as
!              "ghost," but I gotta admit that it's a better way of putting it.
!
! Input:
!   cpGas       Specific heat at constant pressure
!   mmGas       Molecular mass
!   nx,ny,nz    Components of unit normal vector
!   rl          Density
!   ul         x-velocity component
!   vl         y-velocity component
!   wl         z-velocity component
!   fs          Grid speed
!   pl          Pressure
!
! Output: 
!   pl          Pressure
!
! Notes: 
!   1. Valid only for thermally and calorically perfect gas.
!
! ******************************************************************************

      SUBROUTINE RFLU_SetRindStateSlipWallPerf(cpGas,mmGas,nx,ny,nz,rl,
     >                                         ul,vl,wl,fs,pl)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

      real MixtPerf_R_M, MixtPerf_G_CpR, MixtPerf_C_DGP
      external MixtPerf_R_M, MixtPerf_G_CpR, MixtPerf_C_DGP

! ==============================================================================  
!   Arguments 
! ==============================================================================  

      REAL cpGas,fs,mmGas,nx,ny,nz,rl,ul,vl,wl
      REAL pl

! ==============================================================================  
!   Locals 
! ==============================================================================  

      REAL al,gGas,irl,ql,rGas,term
          
! ******************************************************************************
!   Compute wall pressure
! ******************************************************************************

      rGas = MixtPerf_R_M(mmGas)
      gGas = MixtPerf_G_CpR(cpGas,rGas)
 
      irl = 1.0/rl
      ql  = ul*nx + vl*ny + wl*nz - fs
 
      al  = MixtPerf_C_DGP(rl,gGas,pl)

      IF ( ql < 0.0 ) THEN
         term = 1.0 + 0.5*(gGas-1.0)*ql/al
         pl   = pl*term**(2.0*gGas/(gGas-1.0))
      ELSE
         term = (gGas+1.0)/4.0
         pl   = pl + term*rl*ql*(ql + SQRT(al*al+term*term*ql*ql)/term)
      END IF ! ql
 
! ******************************************************************************
!   End
! ******************************************************************************

      END SUBROUTINE RFLU_SetRindStateSlipWallPerf
