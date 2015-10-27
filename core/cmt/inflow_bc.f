c--------------------------------------------------------------------
      subroutine inflow(nvar,f,e,faceq,bcq,flux)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTBCDATA'
      integer nvar,f,e
      real faceq(nvar,nx1,nz1,2*ldim,nelt) 
      real bcq(nvar,nx1,nz1,2*ldim,nelt) 
      real flux(nx1*nz1,2*ldim,nelt,*)  

      ! Assume that the user knows inflow is sub or super sonic

      call inflow_rflu(nvar,f,e,faceq,bcq,flux)
      return
      end
!--------------------------------------------------------------------

      subroutine inflow_rflu(nvar,f,e,faceq,bcq,flux1)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'GEOM'
      include 'PARALLEL'

      integer f,e ! intent(in)
      integer i,bcOptType
      real faceq(nvar,nx1*nz1,2*ldim,nelt) ! intent(in)
      real bcq(nvar,nx1*nz1,2*ldim,nelt)   ! intent(out)
      real flux1(nx1*nz1,2*ldim,nelt,*)    ! intent(out)
      real nx,ny,nz,rl,ul,vl,wl,pl,rr,ur,vr,wr,prf,flx(5),fs
      real ptot,ttot,mach

      nxz=nx1*nz1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg)
         bcOptType=0
         l=l+1
         nx = unx(l,1,f,e)
         ny = uny(l,1,f,e)
         nz = unz(l,1,f,e)
         phl= faceq(iph,l,f,e)
         rl = faceq(iu1,l,f,e)/phl
         rul= faceq(iu2,l,f,e)/phl
         rvl= faceq(iu3,l,f,e)/phl
         rwl= faceq(iu4,l,f,e)/phl
         mach = sqrt(ux**2+uy**2+uz**2)/asnd
         if (mach.lt.1.0) bcOptType=1
         betah = atan2(uy,ux)
         betav = atan2(uz,ux)
! belongs in userbc somehow
         call BcondInflowPerf(bcOptType,0,p0in,t0in,betah,
     >                        betav,mach,nx,ny,nz,cp,molarmass,rl,rul,
     >                        rvl,rwl,rr,rur,rvr,rwr,rer,pres)
         bcq(irho,l,f,e) = rr     ! lol aliased
!        bcq(iux, l,f,e) = rur/rr ! lol aliased
!        bcq(iuy, l,f,e) = rvr/rr ! lol aliased
!        bcq(iuz, l,f,e) = rwr/rr ! lol aliased
         bcq(iux, l,f,e) = ux
         bcq(iuy, l,f,e) = uy
         bcq(iuz, l,f,e) = uz
         bcq(ipr, l,f,e) = pres ! BcondInflowPerf
         bcq(ithm,l,f,e) = temp ! userbc
         bcq(isnd,l,f,e) = asnd ! userbc
         bcq(iph, l,f,e) = phi  ! userbc
         bcq(icvf,l,f,e) = cv   ! userbc
         bcq(icpf,l,f,e) = cp   ! userbc
         bcq(iu1, l,f,e) = phi*rr  ! still aliased
         bcq(iu2, l,f,e) = phi*rur ! still aliased
         bcq(iu3, l,f,e) = phi*rvr ! still aliased
         bcq(iu4, l,f,e) = phi*rwr ! still aliased
         bcq(iu5, l,f,e) = phi*rer ! still aliased
!-----------------------------------------------------------------------
! JH031615 OK fine we need more face storage for gas props, and
!          compute_gas_props_face needs output arguments of some kind
         fs = 0.0 ! moving grid stuff later or never
! JH030915 rewrite the above to use compute_gas_props_face
!-----------------------------------------------------------------------
! Inviscid flux at inflow can probably just be hardcoded instead of
! derived from a trivial call of CentralInviscid_FluxFunction
         pl=faceq(iph,l,f,e)*faceq(ipr,l,f,e) ! needs phi. U has phi
         pres=phi*pres ! needs phi. U in bcq has phi already
         call CentralInviscid_FluxFunction(nx,ny,nz,fs,
     >                                     faceq(iu1,l,f,e),pl,! U-
     >                                     bcq(iu1,l,f,e),pres,! U+,bc
     >                                     flx)
         do j=1,5
            flux1(l,f,e,j)=flx(j) ! this one has phi in it already
         enddo

      enddo
      enddo
      enddo

      return
      end subroutine inflow_rflu

!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: set inflow boundary condition for one cell.
!
! Description: The subsonic boundary condition is based on the extrapolation of 
!   the Riemann invariant from the interior field (Holmes, D.G.: Inviscid 2D 
!   Solutions on Unstructured, Adaptive Grids. VKI Lecture Series 1989-06, 1989). 
!   The supersonic inflow boundary condition computes the conservative variables
!   from given velocity components, density and pressure.
!
! Input: bcOptType  = boundary treatment: subsonic, supersonic, or mixed
!        bcOptFixed = whether _computed_ inflow angle should be fixed or not
!        ptot       = given total pressure
!        ttot       = given total temperature
!        betah      = given inlet angle wrp. to y-axis
!        betav      = given inlet angle wrp. to z-axis
!        sx/y/zn    = components of normalized face vector (outward facing)
!        cpgas      = specific heat at constant pressure (boundary cell)
!        mm         = molecular mass at boundary cell
!        rl         = given density
!        ru/v/wl    = given velocity components
!
! Output: rr      = density at boundary
!         ru/v/wr = density * velocity components at boundary
!         rer     = density * total internal energy at boundary
!         pr      = pressure at boundary
!
! Notes: 
!   1. This condition is valid only for thermally and calorically perfect  
!      gas.
!   2. Important to avoid division by MakeNonZero(sl) when computing eta 
!      because that computation can become undefined for reservoir inflow
!      conditions, i.e., for a vanishing velocity vector. 
!
!******************************************************************************
!
! $Id: BcondInflowPerf.F90,v 1.1.1.1 2014/05/05 21:47:47 tmish Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!******************************************************************************

      SUBROUTINE BcondInflowPerf(bcOptType,bcOptFixed,ptot,ttot,betah,
     >                           betav,mach,sxn,syn,szn,cpgas,mm,rl,rul,
     >                           rvl,rwl,rr,rur,rvr,rwr,rer,pr)
      IMPLICIT NONE

      integer BCOPT_SUBSONIC, BCOPT_MIXED, BCOPT_FIXED_NO
! faked flags
      parameter (BCOPT_SUBSONIC = 1)
      parameter (BCOPT_MIXED = 2)
      parameter (BCOPT_FIXED_NO = 2)
      real  MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT,
     >  MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT,
     >  MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR,
     >  MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo,
     >  MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR,
     >  MixtPerf_T_GMaTo,  MixtPerf_Vm_C2Co2G
      external  MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT,
     >  MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT,
     >  MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR,
     >  MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo,
     >  MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR,
     >  MixtPerf_T_GMaTo,  MixtPerf_Vm_C2Co2G


! ... parameters
      INTEGER bcOptFixed,bcOptType!, INTENT(IN) ::

      REAL betah, betav, cpgas, mach, mm, sxn, syn, szn, !, INTENT(IN) ::
     >                        ptot, rl, rul, rvl, rwl, ttot

      REAL rer, rr, rur, rvr, rwr, pr!, INTENT(OUT) ::

! ... local variables
      REAL al, ar, a02, cp, disc, eta, g, gm1, igm1, ql, rgas, Rm,
     >            sl, sr , tr, ul, ur, vl, vr, wl, wr

!******************************************************************************
! gas properties

      rgas = MixtPerf_R_M(mm)
      g    = MixtPerf_G_CpR(cpgas,rgas)

! subsonic or mixed -----------------------------------------------------------   

      IF ( bcOptType.eq.BCOPT_SUBSONIC.OR.bcOptType.eq.BCOPT_MIXED) THEN
         gm1  = g - 1.0
         igm1 = 1.0/gm1

         ul = rul/rl
         vl = rvl/rl
         wl = rwl/rl

         a02 = MixtPerf_C2_GRT(g,rgas,ttot)

         al = MixtPerf_C_Co2GUVW(a02,g,ul,vl,wl) ! make al consistent with a02
         ql = ul*sxn + vl*syn + wl*szn

! - subsonic

         IF ( ABS(ql) < al ) THEN
            sl = SQRT(ul*ul + vl*vl + wl*wl)
         
            IF ( bcOptFixed .eq. BCOPT_FIXED_NO ) THEN 
               IF ( sl > 1.0E-6 ) THEN ! Avoid ill-defined angle computation
                  eta = ql/sl        
               ELSE 
                  eta = -1.0
               END IF ! sl
            ELSE 
               eta = -1.0
            END IF ! bcOptFixed

            Rm   = al - 0.5*gm1*ql
            disc = 0.5*gm1*eta**2*
     >       (a02/(Rm*Rm)*(1.0 + 0.5*gm1*eta**2) - 1.0)     
        
            IF ( disc < 0.0 ) THEN ! discriminant cannot be negative
               ar = SQRT(a02)
               tr = ttot
               pr = ptot
               sr = 0.0
            ELSE
               ar = Rm/(1.0 + 0.5*gm1*eta*eta)*(1.0+SQRT(disc))
               tr = MixtPerf_T_CGR(ar,g,rgas)
               pr = MixtPerf_P_GPoTTo(g,ptot,tr,ttot)
               sr = MixtPerf_Vm_C2Co2G(ar*ar,a02,g)
            END IF ! disc    
                     
            rr = MixtPerf_D_PRT( pr,rgas,tr )
            ur = sr*COS(betah)*COS(betav)
            vr = sr*SIN(betah)
            wr = sr*COS(betah)*SIN(betav)

            rer = rr*MixtPerf_Eo_DGPVm(rr,g,pr,sr)
            rur = rr*ur
            rvr = rr*vr
            rwr = rr*wr

! - supersonic

         ELSE
            IF ( bcOptType .eq. BCOPT_MIXED ) THEN
               pr = mixtPerf_P_GMaPo(g,mach,ptot)
               tr = mixtPerf_T_GMaTo(g,mach,ttot)
               rr = mixtPerf_D_PRT(pr,rgas,tr)
               ar = mixtPerf_C_GRT(g,rgas,tr)

               ur = mach*ar*COS(betah)*COS(betav)
               vr = mach*ar*SIN(betah)
               wr = mach*ar*COS(betah)*SIN(betav)

               rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
               rur = rr*ur
               rvr = rr*vr
               rwr = rr*wr
            END IF ! bcOptType
         END IF ! ql < al

! supersonic ------------------------------------------------------------------

      ELSE ! bcOptType .eq. BCOPT_SUPERSONIC
         pr = mixtPerf_P_GMaPo(g,mach,ptot)
         tr = mixtPerf_T_GMaTo(g,mach,ttot)
         rr = mixtPerf_D_PRT(pr,rgas,tr)
         ar = mixtPerf_C_GRT(g,rgas,tr)

         ur = mach*ar*COS(betah)*COS(betav)
         vr = mach*ar*SIN(betah)
         wr = mach*ar*COS(betah)*SIN(betav) 

         rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
         rur = rr*ur
         rvr = rr*vr
         rwr = rr*wr
      END IF ! bcOptType

      return
      END SUBROUTINE BcondInflowPerf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondInflowPerf.F90,v $
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
! Revision 1.1  2004/12/01 16:47:56  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/01/29 22:52:36  haselbac
! Added bcOptFixed, fixed bug, clean-up
!
! Revision 1.6  2003/12/04 03:22:56  haselbac
! Fixed bug in formulation, added partial fix for eta
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
