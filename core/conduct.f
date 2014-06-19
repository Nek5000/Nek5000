c-----------------------------------------------------------------------
      subroutine cdscal (igeom)
C
C     Solve the convection-diffusion equation for passive scalar IPSCAL
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      LOGICAL          IFCONV
C
      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT)
     $              ,TB(LX1,LY1,LZ1,LELT)
      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT)
     $              ,H2(LX1,LY1,LZ1,LELT)

      include 'ORTHOT'

      napprox(1) = laxt

      nel    = nelfld(ifield)
      ntot   = nx1*ny1*nz1*nel

      if (igeom.eq.1) then   ! geometry at t^{n-1}

         call makeq
         call lagscal

      else                   ! geometry at t^n

         IF (IFPRINT) THEN
         IF (IFMODEL .AND. IFKEPS) THEN
            NFLDT = NFIELD - 1
            IF (IFIELD.EQ.NFLDT.AND.NID.EQ.0) THEN
               WRITE (6,*) ' Turbulence Model - k/epsilon solution'
            ENDIF
         ELSE
            IF (IFIELD.EQ.2.AND.NID.EQ.0) THEN
               WRITE (6,*) ' Temperature/Passive scalar solution'
            ENDIF
         ENDIF
         ENDIF
         if1=ifield-1
         write(name4,1) if1-1
    1    format('PS',i2)
         if(ifield.eq.2) write(name4,'(A4)') 'TEMP'

C
C        New geometry
C
         isd = 1
         if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2
c        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

         INTYPE = 0
         IF (IFTRAN) INTYPE = -1
         CALL SETHLM  (H1,H2,INTYPE)
         CALL BCNEUSC (TA,-1)
         CALL ADD2    (H2,TA,NTOT)
         CALL BCDIRSC (T(1,1,1,1,IFIELD-1))
         CALL AXHELM  (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMESH,isd)
         CALL SUB3    (TB,BQ(1,1,1,1,IFIELD-1),TA,NTOT)
         CALL BCNEUSC (TA,1)
         CALL ADD2    (TB,TA,NTOT)

c        CALL HMHOLTZ (name4,TA,TB,H1,H2
c    $                 ,TMASK(1,1,1,1,IFIELD-1)
c    $                 ,TMULT(1,1,1,1,IFIELD-1)
c    $                 ,IMESH,TOLHT(IFIELD),NMXH,isd)

         if(iftmsh(ifield)) then
           call hsolve  (name4,TA,TB,H1,H2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approx,napprox,bintm1)
         else
           call hsolve  (name4,TA,TB,H1,H2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approx,napprox,binvm1)
         endif 

         call add2    (t(1,1,1,1,ifield-1),ta,ntot)

         call cvgnlps (ifconv) ! Check convergence for nonlinear problem 
         if (ifconv) goto 2000

C        Radiation case, smooth convergence, avoid flip-flop (ER).
         CALL CMULT (TA,0.5,NTOT)
         CALL SUB2  (T(1,1,1,1,IFIELD-1),TA,NTOT)

 1000    CONTINUE
 2000    CONTINUE
         CALL BCNEUSC (TA,1)
         CALL ADD2 (BQ(1,1,1,1,IFIELD-1),TA,NTOT) ! no idea why... pf

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine makeuq

c     Fill up user defined forcing function and collocate will the
c     mass matrix on the Gauss-Lobatto mesh.

      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'

      ntot = nx1*ny1*nz1*nelfld(ifield)

      time = time-dt        ! Set time to t^n-1 for user function

      call rzero   ( bq(1,1,1,1,ifield-1) ,    ntot)
      call setqvol ( bq(1,1,1,1,ifield-1)          )
      call col2    ( bq(1,1,1,1,ifield-1) ,bm1,ntot)

      time = time+dt        ! Restore time

      return
      end
c-----------------------------------------------------------------------
      subroutine setqvol(bql)

c     Set user specified volumetric forcing function (e.g. heat source).

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      real bql(lx1*ly1*lz1,lelt)

#ifndef MOAB
      nel   = nelfld(ifield)
      nxyz1 = nx1*ny1*nz1
      ntot1 = nxyz1*nel

      do iel=1,nel
         igrp = igroup(iel)
         if (matype(igrp,ifield).eq.1) then ! constant source within a group
            cqvol = cpgrp(igrp,ifield,3)
            call cfill (bql(1,iel),cqvol,nxyz1)
         else  !  pff 2/6/96 ............ default is to look at userq
            call nekuq (bql,iel)
         endif
      enddo
c
c 101 FORMAT(' Wrong material type (',I3,') for group',I3,', field',I2
c    $    ,/,' Aborting in SETQVOL.')
#else
c pulling in temperature right now, since we dont have anything else
      call userq2(bql)
#endif
C   
      return
      end
C
      subroutine nekuq (bql,iel)
C------------------------------------------------------------------
C
C     Generate user-specified volumetric source term (temp./p.s.)
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'PARALLEL'
      include 'TSTEP'
      include 'NEKUSE'
      include 'INPUT'
c
      real bql(lx1,ly1,lz1,lelt)
c
      ielg = lglel(iel)
      do 10 k=1,nz1
      do 10 j=1,ny1
      do 10 i=1,nx1
         call nekasgn (i,j,k,iel)
         qvol = 0.0
         call userq   (i,j,k,ielg)
         bql(i,j,k,iel) = qvol
 10   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine convab
C---------------------------------------------------------------
C
C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C
C---------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
C
      COMMON /SCRUZ/ TA (LX1,LY1,LZ1,LELT)
C
      NEL = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
      CALL CONVOP  (TA,T(1,1,1,1,IFIELD-1))
      CALL COL2    (TA,VTRANS(1,1,1,1,IFIELD),NTOT1)
      CALL SUBCOL3 (BQ(1,1,1,1,IFIELD-1),BM1,TA,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine makeabq
C
C     Sum up contributions to 3rd order Adams-Bashforth scheme.
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRUZ/ TA (LX1,LY1,LZ1,LELT)
C
      AB0   = AB(1)
      AB1   = AB(2)
      AB2   = AB(3)
      NEL   = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
C
      CALL ADD3S2 (TA,VGRADT1(1,1,1,1,IFIELD-1),
     $                VGRADT2(1,1,1,1,IFIELD-1),AB1,AB2,NTOT1)
      CALL COPY   (   VGRADT2(1,1,1,1,IFIELD-1),
     $                VGRADT1(1,1,1,1,IFIELD-1),NTOT1)
      CALL COPY   (   VGRADT1(1,1,1,1,IFIELD-1),
     $                     BQ(1,1,1,1,IFIELD-1),NTOT1)
      CALL ADD2S1 (BQ(1,1,1,1,IFIELD-1),TA,AB0,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine makebdq
C-----------------------------------------------------------------------
C
C     Add contributions to F from lagged BD terms.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
C
      COMMON /SCRNS/ TA (LX1,LY1,LZ1,LELT)
     $ ,             TB (LX1,LY1,LZ1,LELT)
     $ ,             H2 (LX1,LY1,LZ1,LELT)
C
      NEL   = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
      CONST = 1./DT
      CALL COPY  (H2,VTRANS(1,1,1,1,IFIELD),NTOT1)
      CALL CMULT (H2,CONST,NTOT1)
C
      CALL COL3  (TB,BM1,T(1,1,1,1,IFIELD-1),NTOT1)
      CALL CMULT (TB,BD(2),NTOT1)
C
      DO 100 ILAG=2,NBD
         IF (IFGEOM) THEN
            CALL COL3 (TA,BM1LAG(1,1,1,1,ILAG-1),
     $                    TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
         ELSE
            CALL COL3 (TA,BM1,
     $                    TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
         ENDIF
         CALL CMULT (TA,BD(ILAG+1),NTOT1)
         CALL ADD2  (TB,TA,NTOT1)
 100  CONTINUE
C
      CALL COL2 (TB,H2,NTOT1)
      CALL ADD2 (BQ(1,1,1,1,IFIELD-1),TB,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine convch_old
C
C     Compute convective contribution using 
C     operator-integrator-factor method (characteristics).
C
      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRUZ/ TCH (LX1,LY1,LZ1,LELT)
     $ ,             H2  (LX1,LY1,LZ1,LELT)
      LOGICAL IFCHAB
C
      IFCHAB = .FALSE.
      NEL    = NELFLD(IFIELD)
      NTOT1  = NX1*NY1*NZ1*NEL
      CONST  = 1./DT
      CALL COPY  (H2,VTRANS(1,1,1,1,IFIELD),NTOT1)
      CALL CMULT (H2,CONST,NTOT1)
C
      DO 100 ILAG=NBD,1,-1
         IF (IFCHAB) THEN
            CALL THYPAB (TCH,ILAG)
         ELSE
            CALL THYPRK (TCH,ILAG)
         ENDIF
         CALL COL2  (TCH,BM1,NTOT1)
         CALL COL2  (TCH,H2 ,NTOT1)
         CALL CMULT (TCH,BD(ILAG+1),NTOT1)
         CALL ADD2  (BQ(1,1,1,1,IFIELD-1),TCH,NTOT1)
 100  CONTINUE
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine thyprk (tch,ilag)
C
C     Convection of a passive scalar.
C     Runge-Kutta scheme.
C
      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      REAL           TCH   (LX1,LY1,LZ1,1)
      COMMON /SCRNS/ VXN   (LX1,LY1,LZ1,LELV)
     $ ,             VYN   (LX1,LY1,LZ1,LELV)
     $ ,             VZN   (LX1,LY1,LZ1,LELV)
     $ ,             HTMASK(LX1,LY1,LZ1,LELT)
     $ ,             WORK  (LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ RK1   (LX1,LY1,LZ1,LELT)
     $ ,             RK2   (LX1,LY1,LZ1,LELT)
     $ ,             RK3   (LX1,LY1,LZ1,LELT)
     $ ,             RK4   (LX1,LY1,LZ1,LELT)
C
      NTOT1 = NX1*NY1*NZ1*NELV
C
      CALL OPCOPY  (VXN,VYN,VZN,VX,VY,VZ)
      CALL HYPMSK1 (HTMASK)
      CALL TAUINIT (TAU,ILAG)
      CALL TCHINIT (TCH,ILAG)
      CALL VELCONV (VXN,VYN,VZN,TAU)
C
      DO 1000 JLAG=ILAG,1,-1
C
         DTAU   = DTLAG(JLAG)/(NTAUBD)
         DTHALF = DTAU/2.
         CRK1   = DTAU/6.
         CRK2   = DTAU/3.
C
         DO 900 ITAU=1,NTAUBD
C
            CALL FRKCONV (RK1,TCH,HTMASK)
C
            TAU = TAU + DTHALF
            CALL VELCONV (VXN,VYN,VZN,TAU)
C
            CALL COPY    (WORK,TCH,NTOT1)
            CALL ADD2S2  (WORK,RK1,-DTHALF,NTOT1)
            CALL FRKCONV (RK2,WORK,HTMASK)
C
            CALL COPY    (WORK,TCH,NTOT1)
            CALL ADD2S2  (WORK,RK2,-DTHALF,NTOT1)
            CALL FRKCONV (RK3,WORK,HTMASK)
C
            TAU = TAU + DTHALF
            CALL VELCONV (VXN,VYN,VZN,TAU)
C
            CALL COPY    (WORK,TCH,NTOT1)
            CALL ADD2S2  (WORK,RK3,-DTAU,NTOT1)
            CALL FRKCONV (RK4,WORK,HTMASK)
C
            CALL ADD2S2  (TCH,RK1,-CRK1,NTOT1)
            CALL ADD2S2  (TCH,RK2,-CRK2,NTOT1)
            CALL ADD2S2  (TCH,RK3,-CRK2,NTOT1)
            CALL ADD2S2  (TCH,RK4,-CRK1,NTOT1)
C
  900    CONTINUE
 1000 CONTINUE
C
      CALL OPCOPY (VX,VY,VZ,VXN,VYN,VZN)
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine thypab (tch,ilag)
C
C     Convection of a passive scalar.
C     Adams-Bashforth.
C
      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      REAL           TCH   (LX1,LY1,LZ1,1)
      COMMON /SCRNS/ VXN   (LX1,LY1,LZ1,LELV)
     $ ,             VYN   (LX1,LY1,LZ1,LELV)
     $ ,             VZN   (LX1,LY1,LZ1,LELV)
     $ ,             HTMASK(LX1,LY1,LZ1,LELT)
     $ ,             WORK  (LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ TMP1  (LX1,LY1,LZ1,LELT)
     $ ,             TMP2  (LX1,LY1,LZ1,LELT)
C
      NEL   = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
C
      CALL OPCOPY  (VXN,VYN,VZN,VX,VY,VZ)
      CALL HYPMSK1 (HTMASK)
      CALL TAUINIT (TAU,ILAG)
      CALL TCHINIT (TCH,ILAG)
      CALL VELCONV (VXN,VYN,VZN,TAU)
C
      DT2   = 0.
      DT1   = 0.
C
      DO 1000 JLAG=ILAG,1,-1
C
         DT2 = DT1
         DT1 = DTAU
         DTAU   = DTLAG(JLAG)/(NTAUBD)
C
         IF (JLAG.EQ.ILAG) THEN
            DTAU1 = 0.1*DTAU**1.5
            DTAU3 = (DTLAG(JLAG)-2.*DTAU1)/(NTAUBD-2)
         ELSE
            DTAU1 = DTAU 
            DTAU3 = DTAU
         ENDIF
C
         DO 900 ITAU=1,NTAUBD
C
         IF (ITAU.GT.1) THEN
            DT2 = DT1
            DT1 = DTAU
         ENDIF
         IF ((JLAG.EQ.ILAG).AND.(ITAU.LE.2)) THEN
            AB1  = 1.
            AB2  = 0.
            AB3  = 0.
            DTAU = DTAU1
         ELSE
            DTAU = DTAU3
            AB3  = (DTAU/DT2)*((DTAU/3.+DT1/2.)/(DT1+DT2))
            AB2  = -(DTAU/DT1)*(0.5+(DTAU/3.+DT1/2.)/DT2)
            AB1  = 1.-AB2-AB3 
         ENDIF
         TAU = TAU + DTAU
C
         CALL ADD3S2  (WORK,TMP1,TMP2,AB2,AB3,NTOT1)
         CALL COPY    (TMP2,TMP1,NTOT1)
         CALL FRKCONV (TMP1,TCH,HTMASK)
         CALL ADD2S2  (WORK,TMP1,AB1,NTOT1)
         CALL ADD2S2  (TCH,WORK,-DTAU,NTOT1)
C
         CALL VELCONV (VXN,VYN,VZN,TAU)
C
  900    CONTINUE
 1000 CONTINUE
C
      CALL OPCOPY (VX,VY,VZ,VXN,VYN,VZN)
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine hypmsk1 (htmask)
C
C     Generate mask-array for the hyperbolic system (passive scalar).
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'TSTEP'
      REAL           HTMASK(LX1,LY1,LZ1,1)
      CHARACTER      CB*3
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      COMMON /CTMP1/ WORK   (LXYZ1,LELT)
C
      NFACES= 2*NDIM
      NEL   = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
      CALL RZERO(WORK  ,NTOT1)
      CALL RONE (HTMASK,NTOT1)
C
      DO 100 IE=1,NELV
      DO 100 IFACE=1,NFACES
        CB=CBC(IFACE,IE,IFIELD)
        IF (CB.EQ.'T  ' .OR. CB.EQ.'t  ' .OR.
     $      CB.EQ.'KD ' .OR. CB.EQ.'kd ' .OR.
     $      CB.EQ.'ED ' .OR. CB.EQ.'ed ') THEN
         CALL FACCL3 (WORK(1,IE),VX(1,1,1,IE),UNX(1,1,IFACE,IE),IFACE)
         CALL FADDCL3(WORK(1,IE),VY(1,1,1,IE),UNY(1,1,IFACE,IE),IFACE)
         IF (IF3D) 
     $   CALL FADDCL3(WORK(1,IE),VZ(1,1,1,IE),UNZ(1,1,IFACE,IE),IFACE)
         CALL FCAVER (TAVER,WORK,IE,IFACE)
C
         IF (TAVER.LT.0.) CALL FACEV (HTMASK,IE,IFACE,0.0,NX1,NY1,NZ1)
         IF (CB.EQ.'KWS' .OR. CB.EQ.'EWS')
     $   CALL FACEV (HTMASK,IE,IFACE,0.0,NX1,NY1,NZ1)
       ENDIF
 100  CONTINUE
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine tchinit (tch,ilag)
C
C     Set initial conditions for subintegration
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      REAL TCH (LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NELFLD(IFIELD)
      IF (ILAG.EQ.1) THEN
         CALL COPY (TCH,T(1,1,1,1,IFIELD-1),NTOT1)
      ELSE
         CALL COPY (TCH,TLAG(1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
      ENDIF
      return
      end
C
c-----------------------------------------------------------------------
      subroutine lagscal
C-----------------------------------------------------------------------
C
C     Keep old passive scalar field(s) 
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
C
      NTOT1 = NX1*NY1*NZ1*NELFLD(IFIELD)
C
      DO 100 ILAG=NBDINP-1,2,-1
         CALL COPY (TLAG(1,1,1,1,ILAG  ,IFIELD-1),
     $              TLAG(1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
 100  CONTINUE
C
      CALL COPY (TLAG(1,1,1,1,1,IFIELD-1),T(1,1,1,1,IFIELD-1),NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldrq (x,txt10,ichk)
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(nx1,ny1,nz1,lelt)
      character*10 txt10
c
      integer idum,e
      save idum
      data idum /-3/

      if (idum.lt.0) return
c
C
      mtot = nx1*ny1*nz1*nelv
      if (nx1.gt.8.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      nell = nelt
      rnel = nell
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nell-ne+1
      k = 1
      do ie=1,1
         ne = 0
         write(6,116) txt10,k,ie,xmin,xmax,istep,time
         do l=0,1
            write(6,117) 
            do j=ny1,1,-1
              if (nx1.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
              if (nx1.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
              if (nx1.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
              if (nx1.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
              if (nx1.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
              if (nx1.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
              if (nx1.eq.8) write(6,118) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
            enddo
         enddo
      enddo
C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
  118 FORMAT(8f12.9)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
      if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine cdscal_expl (igeom)
C
C     explicit convection-diffusion equation for passive scalar
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      common  /cprint/ ifprint
      logical          ifprint
      logical          ifconv

      common /scrns/ ta(lx1,ly1,lz1,lelt)
     $              ,tb(lx1,ly1,lz1,lelt)
      common /scrvh/ h1(lx1,ly1,lz1,lelt)
     $              ,h2(lx1,ly1,lz1,lelt)


c     QUESTIONABLE support for Robin BC's at this point! (5/15/08)

      nel    = nelfld(ifield)
      ntot   = nx1*ny1*nz1*nel

      if (igeom.eq.1) then   ! geometry at t^{n-1}

         call makeq
         call lagscal

      else                   ! geometry at t^n

         if (.true..and.nio.eq.0) 
     $      write (6,*) istep,ifield,' explicit step'


C        New geometry

         isd = 1
         if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         intype = 0
         if (iftran) intype = -1
         call sethlm  (h1,h2,intype)

         call bcneusc (ta,-1)       ! Modify diagonal for Robin condition
         call add2    (h2,ta ,ntot)
         call col2    (h2,BM1,ntot)

         call bcneusc (tb,1)        ! Modify rhs for flux bc
         call add2    (bq(1,1,1,1,ifield-1),tb,ntot)

         call dssum   (bq(1,1,1,1,ifield-1),nx1,ny1,nz1)
         call dssum   (h2,nx1,ny1,nz1)

         call invcol3 (t(1,1,1,1,ifield-1),bq(1,1,1,1,ifield-1),h2,ntot)

         call bcdirsc (t(1,1,1,1,ifield-1)) ! --> no mask needed

      endif                   ! geometry at t^n

      return
      end
c-----------------------------------------------------------------------
      subroutine diffab  ! explicit treatment of diffusion operator
c
c     Eulerian scheme, add diffusion term to forcing function 
c     at current time step.
c

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'INPUT'

      common /scruz/ ta(lx1,ly1,lz1,lelt)
     $              ,h2(lx1,ly1,lz1,lelt)

      nel = nelfld(ifield)
      ntot1 = nx1*ny1*nz1*nel

      intype = 0
      if (iftran) intype = -1

      isd = 1
      if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

      imesh = 1
c      if (iftmsh(ifield)) imesh=2

      call rzero   (h2,ntot1)
      call axhelm  (ta,t(1,1,1,1,ifield-1),vdiff(1,1,1,1,ifield)
     $             ,h2,imesh,isd)
      call sub2    (bq(1,1,1,1,ifield-1),ta,ntot1)

      return
      end
c-----------------------------------------------------------------------
