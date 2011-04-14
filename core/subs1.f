c-----------------------------------------------------------------------
      SUBROUTINE SETDT
c
c     Set the new time step. All cases covered.
c
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      COMMON /CPRINT/ IFPRINT
      LOGICAL         IFPRINT
      COMMON /UDXMAX/ UMAX
      REAL     DTOLD
      SAVE     DTOLD
      DATA     DTOLD /0.0/
      REAL     DTOpf
      SAVE     DTOpf
      DATA     DTOpf /0.0/
      logical iffxdt
      save    iffxdt
      data    iffxdt /.false./
C
      if (param(12).lt.0.or.iffxdt) then
         iffxdt    = .true.
         param(12) = abs(param(12))
         dt        = param(12)
         dtopf     = dt 
         call compute_cfl(umax,vx,vy,vz,1.0)
         goto 200
      else IF (PARAM(84).NE.0.0) THEN
         if (dtold.eq.0.0) then
            dt   =param(84)
            dtold=param(84)
            dtopf=param(84)
            return
         else
            dtold=dt
            dtopf=dt
            dt=dtopf*param(85)
            dt=min(dt,param(12))
         endif
      ENDIF
C
C     Find DT=DTCFL based on CFL-condition (if applicable)
C
      CALL SETDTC
      DTCFL = DT
C
C     Find DTFS based on surface tension (if applicable)
C
      CALL SETDTFS (DTFS)
C
C     Select appropriate DT
C
      IF ((DT.EQ.0.).AND.(DTFS.GT.0.)) THEN
          DT = DTFS
      ELSEIF ((DT.GT.0.).AND.(DTFS.GT.0.)) THEN
          DT = MIN(DT,DTFS)
      ELSEIF ((DT.EQ.0.).AND.(DTFS.EQ.0.)) THEN
          DT = 0.
          IF (IFFLOW.AND.NID.EQ.0.AND.IFPRINT) THEN
             WRITE (6,*) 'WARNING: CFL-condition & surface tension'
             WRITE (6,*) '         are not applicable'
          ENDIF
      ELSEIF ((DT.GT.0.).AND.(DTFS.EQ.0.)) THEN
          DT = DT
      ELSE
          DT = 0.
          IF (NID.EQ.0) WRITE (6,*) 'WARNING: DT<0 or DTFS<0'
          IF (NID.EQ.0) WRITE (6,*) '         Reset DT      '
      ENDIF
C
C     Check DT against user-specified input, DTINIT=PARAM(12).
C
      IF ((DT.GT.0.).AND.(DTINIT.GT.0.)) THEN
         DT = MIN(DT,DTINIT)
      ELSEIF ((DT.EQ.0.).AND.(DTINIT.GT.0.)) THEN
         DT = DTINIT
      ELSEIF ((DT.GT.0.).AND.(DTINIT.EQ.0.)) THEN
         DT = DT
      ELSEIF (.not.iffxdt) THEN
         DT = 0.001
         IF(NID.EQ.0)WRITE (6,*) 'WARNING: Set DT=0.001 (arbitrarily)'
      ENDIF
C
C     Check if final time (user specified) has been reached.
C
 200  IF (TIME+DT .GE. FINTIM .AND. FINTIM.NE.0.0) THEN
C        Last step
         LASTEP = 1
         DT = FINTIM-TIME
         IF (NID.EQ.0) WRITE (6,*) 'Final time step = ',DT
      ENDIF
C
      COURNO = DT*UMAX
      IF (NID.EQ.0.AND.IFPRINT.AND.DT.NE.DTOLD)
     $   WRITE (6,100) DT,DTCFL,DTFS,DTINIT
 100     FORMAT(5X,'DT/DTCFL/DTFS/DTINIT',4E12.3)
C
C     Put limits on how much DT can change.
C
      IF (DTOLD.NE.0.0) THEN
         DTMIN=0.8*DTOLD
         DTMAX=1.2*DTOLD
         DT = MIN(DTMAX,DT)
         DT = MAX(DTMIN,DT)
      ENDIF
      DTOLD=DT

C      IF (PARAM(84).NE.0.0) THEN
C            dt=dtopf*param(85)
C            dt=min(dt,param(12))
C      ENDIF

      if (iffxdt) dt=dtopf
      COURNO = DT*UMAX
c
      if (iffxdt.and.abs(courno).gt.10.*abs(ctarg)) then
         if (nid.eq.0) write(6,*) 'CFL, Ctarg!',courno,ctarg
         call emerxit
      endif

c     Synchronize time step for multiple sessions
      if (IFNEKNEK) then 
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
         DT=glmin(DT,1)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original
       end if  


      RETURN
      END
C
C--------------------------------------------------------
C
      SUBROUTINE CVGNLPS (IFCONV)
C----------------------------------------------------------------------
C
C     Check convergence for non-linear passisve scalar solver.
C     Relevant for solving heat transport problems with radiation b.c.
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      LOGICAL  IFCONV
C
      IF (IFNONL(IFIELD)) THEN
         IFCONV = .FALSE.
      ELSE
         IFCONV = .TRUE.
         RETURN
      ENDIF
C
      TNORM1 = TNRMH1(IFIELD-1)
      CALL UNORM
      TNORM2 = TNRMH1(IFIELD-1)
      EPS = ABS((TNORM2-TNORM1)/TNORM2)
      IF (EPS .LT. TOLNL) IFCONV = .TRUE.
C
      RETURN
      END
C
      SUBROUTINE UNORM
C---------------------------------------------------------------------
C
C     Norm calculation.
C
C--------------------------------------------------------------------- 
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      IF (IFIELD.EQ.1) THEN
C
C        Compute norms of the velocity.
C        Compute time mean (L2) of the inverse of the time step. 
C        Compute L2 in time, H1 in space of the velocity.
C
         CALL NORMVC (VNRMH1,VNRMSM,VNRML2,VNRML8,VX,VY,VZ)
         IF (ISTEP.EQ.0) RETURN
         IF (ISTEP.EQ.1) THEN
            DTINVM = 1./DT
            VMEAN  = VNRML8
         ELSE
            tden   = time
            if (time.le.0) tden = abs(time)+1.e-9
            arg    = ((TIME-DT)*DTINVM**2+1./DT)/tden
            if (arg.gt.0) DTINVM = SQRT(arg)
            arg    = ((TIME-DT)*VMEAN**2+DT*VNRMH1**2)/tden
            if (arg.gt.0) VMEAN  = SQRT(arg)
         ENDIF
      ELSE
C
C     Compute norms of a passive scalar
C
         CALL NORMSC (TNRMH1(IFIELD-1),TNRMSM(IFIELD-1),
     $                TNRML2(IFIELD-1),TNRML8(IFIELD-1),
     $                     T(1,1,1,1,IFIELD-1),IMESH)
         TMEAN(IFIELD-1) = 0.
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE CHKTMG (TOL,RES,W1,W2,MULT,MASK,IMESH)
C-------------------------------------------------------------------
C
C     Check that the tolerances are not too small for the MG-solver.
C     Important when calling the MG-solver (Gauss-Lobatto mesh).
C     Note: direct stiffness summation
C
C-------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'EIGEN'
C
      REAL RES  (LX1,LY1,LZ1,1)
      REAL W1   (LX1,LY1,LZ1,1)
      REAL W2   (LX1,LY1,LZ1,1)
      REAL MULT (LX1,LY1,LZ1,1)
      REAL MASK (LX1,LY1,LZ1,1)
C
C     Single or double precision???
C
      DELTA = 1.E-9
      X     = 1.+DELTA
      Y     = 1.
      DIFF  = ABS(X-Y)
      IF (DIFF.EQ.0.) EPS = 1.E-6*EIGGA/EIGAA
      IF (DIFF.GT.0.) EPS = 1.E-13*EIGGA/EIGAA
C
      IF (IMESH.EQ.1) NL = NELV
      IF (IMESH.EQ.2) NL = NELT
      NTOT1 = NX1*NY1*NZ1*NL
      CALL COPY (W1,RES,NTOT1)
C
      CALL DSSUM (W1,NX1,NY1,NZ1)
C
      IF (IMESH.EQ.1) THEN
         CALL COL3 (W2,BINVM1,W1,NTOT1)
         RINIT  = SQRT(GLSC3 (W2,W1,MULT,NTOT1)/VOLVM1)
      ELSE
         CALL COL3 (W2,BINTM1,W1,NTOT1)
         RINIT  = SQRT(GLSC3 (W2,W1,MULT,NTOT1)/VOLTM1)
      ENDIF
      RMIN   = EPS*RINIT
      IF (TOL.LT.RMIN) THEN
         TOLOLD=TOL
         TOL = RMIN
         IF (NID.EQ.0)
     $   WRITE (6,*) 'New MG-tolerance (RINIT*epsm*cond) = ',TOL,TOLOLD
      ENDIF
C
      CALL RONE (W1,NTOT1)
      BCNEU1 = GLSC3(W1,MASK,MULT,NTOT1)
      BCNEU2 = GLSC3(W1,W1  ,MULT,NTOT1)
      BCTEST = ABS(BCNEU1-BCNEU2)
      IF (BCTEST .LT. .1) THEN
         OTR = GLSUM (RES,NTOT1)
         TOLMIN = ABS(OTR)*EIGGA/EIGAA
         IF (TOL .LT. TOLMIN) THEN
             TOLOLD = TOL
             TOL = TOLMIN
             IF (NID.EQ.0)
     $       WRITE (6,*) 'New MG-tolerance (OTR) = ',TOL,TOLOLD
         ENDIF
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE SETDTC
C--------------------------------------------------------------
C
C     Compute new timestep based on CFL-condition
C
C--------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'MVGEOM'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      COMMON /CTMP1/ U(LX1,LY1,LZ1,LELV)
     $ ,             V(LX1,LY1,LZ1,LELV)
     $ ,             W(LX1,LY1,LZ1,LELV)
      COMMON /CTMP0/ X(LX1,LY1,LZ1,LELV)
     $ ,             R(LX1,LY1,LZ1,LELV)
      COMMON /UDXMAX/ UMAX
C
      REAL WORK(1)
C
      REAL VCOUR
      SAVE VCOUR
C
      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST/0/
C
C
C     Steady state => all done
C
C
      IF (.NOT.IFTRAN) THEN
         IFIRST=1
         LASTEP=1
         RETURN
      ENDIF
C
      irst = param(46)
      if (irst.gt.0) ifirst=1
C
C     First time around
C
      IF (IFIRST.EQ.0) THEN
         DT     = DTINIT
         IF (IFFLOW) THEN
            IFIELD = 1
            CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
         ENDIF
      ENDIF
      IFIRST=IFIRST+1
C
      DTOLD = DT
C
C     Convection ?
C
C     Don't enforce Courant condition if there is no convection.
C
C
      ICONV=0
      IF (IFFLOW .AND. IFNAV) ICONV=1
      IF (IFWCNO)             ICONV=1
      IF (IFHEAT) THEN
         DO 10 IPSCAL=0,NPSCAL
            IF (IFADVC(IPSCAL+2)) ICONV=1
   10    CONTINUE
      ENDIF
      IF (ICONV.EQ.0) THEN
         DT=0.
         RETURN
      ENDIF
C
C
C     Find Courant and Umax
C
C
      NTOT   = NX1*NY1*NZ1*NELV
      NTOTL  = LX1*LY1*LZ1*LELV
      NTOTD  = NTOTL*NDIM
      COLD   = COURNO
      CMAX   = 1.2*CTARG
      CMIN   = 0.8*CTARG
      CALL CUMAX (VX,VY,VZ,UMAX)
C
C     Zero DT
C
      IF (DT .EQ. 0.0) THEN
C
         IF (UMAX .NE. 0.0) THEN
            DT = CTARG/UMAX
            VCOUR = UMAX
         ELSEIF (IFFLOW) THEN
C
C           We'll use the body force to predict max velocity
C
            CALL SETPROP
            IFIELD = 1
C
            CALL MAKEUF
            CALL OPDSSUM (BFX,BFY,BFZ)
            CALL OPCOLV  (BFX,BFY,BFZ,BINVM1)
            FMAX=0.0
            CALL RZERO (U,NTOTD)
            DO 600 I=1,NTOT
               U(I,1,1,1) = ABS(BFX(I,1,1,1))
               V(I,1,1,1) = ABS(BFY(I,1,1,1))
               W(I,1,1,1) = ABS(BFZ(I,1,1,1))
 600        CONTINUE
            FMAX    = GLMAX (U,NTOTD)
            DENSITY = AVTRAN(1)
            AMAX    = FMAX/DENSITY
            WORK(1) = SQRT( (XM1(1,1,1,1)-XM1(2,1,1,1))**2 +
     $                      (YM1(1,1,1,1)-YM1(2,1,1,1))**2 + 
     $                      (ZM1(1,1,1,1)-ZM1(2,1,1,1))**2 ) 
            DXCHAR  = GLMIN (WORK,1)
            IF (AMAX.NE.0.) THEN
               DT = SQRT(CTARG*DXCHAR/AMAX)
            ELSE
               IF (NID.EQ.0) 
     $         WRITE (6,*) 'CFL: Zero velocity and body force'
               DT = 0.0
               RETURN
            ENDIF
         ELSEIF (IFWCNO) THEN
            IF (NID.EQ.0) 
     $      WRITE (6,*) ' Stefan problem with no fluid flow'
            DT = 0.0
            RETURN
         ENDIF
C
      ELSEIF ((DT.GT.0.0).AND.(UMAX.NE.0.0)) THEN
C
C
C     Nonzero DT & nonzero velocity
C
C
      COURNO = DT*UMAX
      VOLD   = VCOUR
      VCOUR  = UMAX
      IF (IFIRST.EQ.1) THEN
         COLD = COURNO
         VOLD = VCOUR
      ENDIF
      CPRED  = 2.*COURNO-COLD
C
C     Change DT if it is too big or if it is too small 
C
c     if (nid.eq.0) 
c    $write(6,917) dt,umax,vold,vcour,cpred,cmax,courno,cmin
c 917 format(' dt',4f9.5,4f10.6)
      IF(COURNO.GT.CMAX .OR. CPRED.GT.CMAX .OR. COURNO.LT.CMIN) THEN
C
            A=(VCOUR-VOLD)/DT
            B=VCOUR
C           -C IS Target Courant number
            C=-CTARG
            DISCR=B**2-4*A*C
            DTOLD=DT
            IF(DISCR.LE.0.0)THEN
               if (nid.eq.0) 
     $         PRINT*,'Problem calculating new DT Discriminant=',discr
               DT=DT*(CTARG/COURNO)
C               IF(DT.GT.DTOLD) DT=DTOLD
            ELSE IF(ABS((VCOUR-VOLD)/VCOUR).LT.0.001)THEN
C              Easy: same v as before (LINEARIZED)
               DT=DT*(CTARG/COURNO)
c     if (nid.eq.0) 
c    $write(6,918) dt,dthi,dtlow,discr,a,b,c
c 918 format(' d2',4f9.5,4f10.6)
            ELSE
               DTLOW=(-B+SQRT(DISCR) )/(2.0*A)
               DTHI =(-B-SQRT(DISCR) )/(2.0*A)
               IF(DTHI .GT. 0.0 .AND. DTLOW .GT. 0.0)THEN
                  DT = MIN (DTHI,DTLOW)
c     if (nid.eq.0) 
c    $write(6,919) dt,dthi,dtlow,discr,a,b,c
c 919 format(' d3',4f9.5,4f10.6)
               ELSE IF(DTHI .LE. 0.0 .AND. DTLOW .LE. 0.0)THEN
c                 PRINT*,'DTLOW,DTHI',DTLOW,DTHI
c                 PRINT*,'WARNING: Abnormal DT from CFL-condition'
c                 PRINT*,'         Keep going'
                  DT=DT*(CTARG/COURNO)
               ELSE
C                 Normal case; 1 positive root, one negative root
                  DT = MAX (DTHI,DTLOW)
c     if (nid.eq.0) 
c    $write(6,929) dt,dthi,dtlow,discr,a,b,c
c 929 format(' d4',4f9.5,4f10.6)
               ENDIF
            ENDIF
C           We'll increase gradually-- make it the geometric mean between
c     if (nid.eq.0) 
c    $write(6,939) dt,dtold
c 939 format(' d5',4f9.5,4f10.6)
            IF (DTOLD/DT .LT. 0.2) DT = DTOLD*5
      ENDIF
C
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE CUMAX (V1,V2,V3,UMAX)
C
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
C
      COMMON /SCRNS/ XRM1 (LX1,LY1,LZ1,LELV)
     $ ,             XSM1 (LX1,LY1,LZ1,LELV)
     $ ,             XTM1 (LX1,LY1,LZ1,LELV)
     $ ,             YRM1 (LX1,LY1,LZ1,LELV)
     $ ,             YSM1 (LX1,LY1,LZ1,LELV)
     $ ,             YTM1 (LX1,LY1,LZ1,LELV)
      COMMON /SCRMG/ ZRM1 (LX1,LY1,LZ1,LELV)
     $ ,             ZSM1 (LX1,LY1,LZ1,LELV)
     $ ,             ZTM1 (LX1,LY1,LZ1,LELV)
      COMMON /CTMP1/ U    (LX1,LY1,LZ1,LELV)
     $ ,             V    (LX1,LY1,LZ1,LELV)
     $ ,             W    (LX1,LY1,LZ1,LELV)
      COMMON /CTMP0/ X    (LX1,LY1,LZ1,LELV)
     $ ,             R    (LX1,LY1,LZ1,LELV)
      COMMON /DELRST/ DRST(LX1),DRSTI(LX1)
C
      DIMENSION V1(LX1,LY1,LZ1,1)
     $        , V2(LX1,LY1,LZ1,1)
     $        , V3(LX1,LY1,LZ1,1)
      DIMENSION U3(3)
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /0/
C
      NTOT  = NX1*NY1*NZ1*NELV
      NTOTL = LX1*LY1*LZ1*LELV
      NTOTD = NTOTL*NDIM
C
C     Compute isoparametric partials.
C
      CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,
     $             IFAXIS)
C
C     Compute maximum U/DX
C
      IF (ICALLD.EQ.0) THEN
         ICALLD=1
         DRST (1)=ABS(ZGM1(2,1)-ZGM1(1,1))
         DRSTI(1)=1.0/DRST(1)
         DO 400 I=2,NX1-1
            DRST (I)=ABS(ZGM1(I+1,1)-ZGM1(I-1,1))/2.0
            DRSTI(I)=1.0/DRST(I)
 400     CONTINUE
         DRST (NX1)=DRST(1)
         DRSTI(NX1)=1.0/DRST(NX1)
      ENDIF
C
C     Zero out scratch arrays U,V,W for ALL declared elements...
C
      CALL RZERO3 (U,V,W,NTOTL)
C
      IF (NDIM.EQ.2) THEN

      CALL VDOT2  (U,V1  ,V2  ,RXM1,RYM1,NTOT)
      CALL VDOT2  (R,RXM1,RYM1,RXM1,RYM1,NTOT)
      CALL VDOT2  (X,XRM1,YRM1,XRM1,YRM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(U,R,NTOT)
C
      CALL VDOT2  (V,V1  ,V2  ,SXM1,SYM1,NTOT)
      CALL VDOT2  (R,SXM1,SYM1,SXM1,SYM1,NTOT)
      CALL VDOT2  (X,XSM1,YSM1,XSM1,YSM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(V,R,NTOT)
C
      ELSE
C
      CALL VDOT3  (U,V1  ,V2  ,V3  ,RXM1,RYM1,RZM1,NTOT)
      CALL VDOT3  (R,RXM1,RYM1,RZM1,RXM1,RYM1,RZM1,NTOT)
      CALL VDOT3  (X,XRM1,YRM1,ZRM1,XRM1,YRM1,ZRM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(U,R,NTOT)
C
      CALL VDOT3  (V,V1  ,V2  ,V3  ,SXM1,SYM1,SZM1,NTOT)
      CALL VDOT3  (R,SXM1,SYM1,SZM1,SXM1,SYM1,SZM1,NTOT)
      CALL VDOT3  (X,XSM1,YSM1,ZSM1,XSM1,YSM1,ZSM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(V,R,NTOT)
C
      CALL VDOT3  (W,V1  ,V2  ,V3  ,TXM1,TYM1,TZM1,NTOT)
      CALL VDOT3  (R,TXM1,TYM1,TZM1,TXM1,TYM1,TZM1,NTOT)
      CALL VDOT3  (X,XTM1,YTM1,ZTM1,XTM1,YTM1,ZTM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(W,R,NTOT)
C
      ENDIF
C
      DO 500 IE=1,NELV
      DO 500 IX=1,NX1
      DO 500 IY=1,NY1
      DO 500 IZ=1,NZ1
            U(IX,IY,IZ,IE)=ABS( U(IX,IY,IZ,IE)*DRSTI(IX) )
            V(IX,IY,IZ,IE)=ABS( V(IX,IY,IZ,IE)*DRSTI(IY) )
            W(IX,IY,IZ,IE)=ABS( W(IX,IY,IZ,IE)*DRSTI(IZ) )
  500    CONTINUE
C
      U3(1)   = VLMAX(U,NTOT)
      U3(2)   = VLMAX(V,NTOT)
      U3(3)   = VLMAX(W,NTOT)
      UMAX    = GLMAX(U3,3)
C
      RETURN
      END
C
      SUBROUTINE SETDTFS (DTFS)
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      COMMON /CTMP0/  STC(LX1,LY1,LZ1),SIGST(LX1,LY1)
      COMMON /SCRVH/  DTST(LX1,LY1)
      CHARACTER CB*3,CB2*2
      REAL WORK(1)
C
C     Applicable?
C
      IF (.NOT.IFSURT) THEN
         DTFS = 0.0
         RETURN
      ENDIF
C
      NFACE  = 2*NDIM
      NXZ1   = NX1*NZ1
      NTOT1  = NX1*NY1*NZ1*NELV
      DTFS   = 1.1E+10
C
C     Kludge : for debugging purpose Fudge Factor comes in
C              as PARAM(45)
C
      FACTOR = PARAM(45)
      IF (FACTOR .LE. 0.0) FACTOR=1.0
      FACTOR = FACTOR / SQRT( PI**3 )
C
      IF (ISTEP.EQ.1) CALL SETPROP
C
      IFIELD = 1
      RHOMIN = GLMIN( VTRANS(1,1,1,1,IFIELD),NTOT1 )
C
      DO 100 IEL=1,NELV
      DO 100 IFC=1,NFACE
         CB =CBC(IFC,IEL,IFIELD)
         CB2=CBC(IFC,IEL,IFIELD)
         IF (CB2.NE.'MS' .AND. CB2.NE.'ms') GOTO 100
             IF (CB2(1:1).EQ.'M') THEN
                BC4 = BC(4,IFC,IEL,IFIELD)
                CALL CFILL  (SIGST,BC4,NXZ1)
             ELSE
                CALL FACEIS (CB,STC,IEL,IFC,NX1,NY1,NZ1)
                CALL FACEXS (SIGST,STC,IFC,0)
             ENDIF
             SIGMAX = VLMAX (SIGST,NXZ1)
             IF (SIGMAX.LE.0.0) GOTO 100
             RHOSIG = SQRT( RHOMIN/SIGMAX )
             IF (NDIM.EQ.2) THEN
                CALL CDXMIN2 (DTST,RHOSIG,IEL,IFC,IFAXIS)
             ELSE
                CALL CDXMIN3 (DTST,RHOSIG,IEL,IFC)
             ENDIF
             DTMIN = VLMIN( DTST,NXZ1 )
             DTFS  = MIN( DTFS,DTMIN )
  100 CONTINUE
      WORK(1) = DTFS
      DTFS    = GLMIN( WORK,1 )
C
      IF (DTFS .GT. 1.E+10) DTFS = 0.0
      DTFS = DTFS * FACTOR
C
      IF (DTFS.EQ.0.0) THEN
         IF (ISTEP.EQ.1.AND.NID.EQ.0) THEN
            WRITE (6,*) ' Warning - zero surface-tension may results in'
            WRITE (6,*) ' instability of free-surface update'
         ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE CDXMIN2 (DTST,RHOSIG,IEL,IFC,IFAXIS)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'DXYZ'
      COMMON /DELRST/ DRST(LX1),DRSTI(LX1)
      COMMON /CTMP0/  XFM1(LX1),YFM1(LX1),T1XF(LX1),T1YF(LX1)
      DIMENSION DTST(LX1,1)
      LOGICAL IFAXIS
C
      DELTA = 1.E-9
      X     = 1.+DELTA
      Y     = 1.
      DIFF  = ABS(X-Y)
      IF (DIFF.EQ.0.) EPS = 1.E-6
      IF (DIFF.GT.0.) EPS = 1.E-13
C
      CALL FACEC2 (XFM1,YFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL),IFC)
C
      IF (IFC.EQ.1 .OR. IFC.EQ.3) THEN
         CALL MXM (DXM1,NX1,XFM1,NX1,T1XF,1)
         CALL MXM (DXM1,NX1,YFM1,NX1,T1YF,1)
      ELSE
         IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
         CALL MXM (DYM1,NY1,XFM1,NY1,T1XF,1)
         CALL MXM (DYM1,NY1,YFM1,NY1,T1YF,1)
      ENDIF
C
      IF (IFAXIS) THEN
         DO 100 IX=1,NX1
            IF (YFM1(IX) .LT. EPS) THEN
               DTST(IX,1) = 1.e+10
            ELSE
               XJ = SQRT( T1XF(IX)**2 + T1YF(IX)**2 )*DRST(IX)
               DTST(IX,1) = RHOSIG * SQRT( XJ**3 ) * YFM1(IX)
            ENDIF
  100    CONTINUE
      ELSE
         DO 200 IX=1,NX1
            XJ = SQRT( T1XF(IX)**2 + T1YF(IX)**2 )*DRST(IX)
            DTST(IX,1) = RHOSIG * SQRT( XJ**3 )
  200    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE CDXMIN3 (DTST,RHOSIG,IEL,IFC)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'DXYZ'
      COMMON /DELRST/ DRST(LX1),DRSTI(LX1)
      COMMON /CTMP0/  XFM1(LX1,LY1),YFM1(LX1,LY1),ZFM1(LX1,LY1)
      COMMON /CTMP1/  DRM1(LX1,LX1),DRTM1(LX1,LY1)
     $             ,  DSM1(LX1,LX1),DSTM1(LX1,LY1)
      COMMON /SCRMG/  XRM1(LX1,LY1),YRM1(LX1,LY1),ZRM1(LX1,LY1)
     $             ,  XSM1(LX1,LY1),YSM1(LX1,LY1),ZSM1(LX1,LY1)
      DIMENSION DTST(LX1,LY1)
C
      CALL FACEXV (XFM1,YFM1,ZFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL),
     $             ZM1(1,1,1,IEL),IFC,0)
      CALL SETDRS (DRM1,DRTM1,DSM1,DSTM1,IFC)
C
      CALL MXM (DRM1,NX1, XFM1,NX1,XRM1,NY1)
      CALL MXM (DRM1,NX1, YFM1,NX1,YRM1,NY1)
      CALL MXM (DRM1,NX1, ZFM1,NX1,ZRM1,NY1)
      CALL MXM (XFM1,NX1,DSTM1,NY1,XSM1,NY1)
      CALL MXM (YFM1,NX1,DSTM1,NY1,YSM1,NY1)
      CALL MXM (ZFM1,NX1,DSTM1,NY1,ZSM1,NY1)
C
      DO 100 IX=1,NX1
      DO 100 IY=1,NY1
         DELR = XRM1(IX,IY)**2 + YRM1(IX,IY)**2 + ZRM1(IX,IY)**2
         DELS = XSM1(IX,IY)**2 + YSM1(IX,IY)**2 + ZSM1(IX,IY)**2
         DELR = SQRT( DELR )*DRST(IX)
         DELS = SQRT( DELS )*DRST(IY)
         XJ   = MIN( DELR,DELS )
         DTST(IX,IY) = RHOSIG * SQRT( XJ**3 )
  100 CONTINUE
C
      RETURN
      END
C
      FUNCTION FACDOT(A,B,IFACE1)
C
C     Take the dot product of A and B on the surface IFACE1 of element IE.
C
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         5 Jan 1989 15:12:22      PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      SUM=0.0
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         SUM = SUM + A(J1,J2,1)*B(I,1)
  100 CONTINUE
C
      FACDOT = SUM
C
      RETURN
      END
C
      SUBROUTINE FCAVER(XAVER,A,IEL,IFACE1)
C------------------------------------------------------------------------
C
C     Compute the average of A over the face IFACE1 in element IEL.
C
C         A is a (NX,NY,NZ) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TOPOL'
      REAL A(LX1,LY1,LZ1,1)
C
      FCAREA = 0.
      XAVER  = 0.
C
C     Set up counters
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         FCAREA = FCAREA+AREA(I,1,IFACE1,IEL)
         XAVER  = XAVER +AREA(I,1,IFACE1,IEL)*A(J1,J2,1,IEL)
  100 CONTINUE
C
      XAVER = XAVER/FCAREA
      RETURN
      END
      SUBROUTINE FACCL2(A,B,IFACE1)
C
C     Collocate B with A on the surface IFACE1 of element IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         5 Jan 1989 15:12:22      PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = A(J1,J2,1)*B(I,1)
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE FACCL3(A,B,C,IFACE1)
C
C     Collocate B with A on the surface IFACE1 of element IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         5 Jan 1989 15:12:22      PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = B(J1,J2,1)*C(I,1)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE FADDCL3(A,B,C,IFACE1)
C
C     Collocate B with C and add to A on the surface IFACE1 of element IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,NZ) data structure
C         C is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         29 Jan 1990 18:00 PST   PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = A(J1,J2,1) + B(J1,J2,1)*C(I,1)
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine sethlm (h1,h2,intloc)
 
c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     (associated with variable IFIELD)
c     INTLOC =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      real h1(1),h2(1)

      nel   = nelfld(ifield)
      ntot1 = nx1*ny1*nz1*nel

      if (iftran) then
         dtbd = bd(1)/dt
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         if (intloc.eq.0) then
            call rzero (h2,ntot1)
         else
            if (ifield.eq.1.or.param(107).eq.0) then 

               call cmult2 (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

            else   ! unsteady reaction-diffusion type equation

               do i=1,ntot1
                 h2(i) = dtbd*vtrans(i,1,1,1,ifield) + param(107)
               enddo

            endif

         endif

c        if (ifield.eq.1 .and. ifanls) then   ! this should be replaced
c           const = 2.                        ! with a correct stress
c           call cmult (h1,const,ntot1)       ! formulation
c        endif

      ELSE
         CALL COPY  (H1,VDIFF (1,1,1,1,IFIELD),NTOT1)
         CALL RZERO (H2,NTOT1)
         if (param(107).ne.0) then
            write(6,*) 'SPECIAL SETHLM!!',param(107)
c           call cfill (h2,param(107),ntot1)
            call copy  (h2,vtrans(1,1,1,1,ifield),ntot1)
         endif
      ENDIF

      return
      end
c-----------------------------------------------------------------------

      SUBROUTINE VPROPS 
C-----------------------------------------------------------------------
C
C     Set material properties
C
C     Material type: 0 for default  (PARAM and PCOND/PRHOCP)
C                    1 for constant props; 
C                    2 for fortran function;
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      LOGICAL  IFKFLD,IFEFLD
C
      NXYZ1 = NX1*NY1*NZ1
      NEL   = NELFLD(IFIELD)
      NTOT1 = NXYZ1*NEL
C
      IF (ISTEP.EQ.0) THEN
C
C        First time around, set defaults
C
         ifvarp(ifield) = .false.
         if (iflomach) ifvarp(ifield) = .true.

         if (.not.ifvarp(ifield)) then ! check all groups
            do iel=1,nel
               igrp  = igroup(iel)
               itype = matype(igrp,ifield)
               if(itype.ne.0) ifvarp(ifield) = .true.
            enddo
         endif

         itest = 0                        ! test against all processors
         if (ifvarp(ifield)) itest = 1
         itest = iglmax(itest,1)
         if (itest.gt.0) ifvarp(ifield) = .true.

      endif         
C
C     Fill up property arrays every time step
C
C     First, check for turbulence models
C
      IF (IFMODEL .AND. IFKEPS) THEN
         CALL TURBFLD (IFKFLD,IFEFLD)
         IF (IFKFLD)           CALL TPROPK
         IF (IFEFLD)           CALL TPROPE
         IF (IFKFLD.OR.IFEFLD) RETURN
      ENDIF
C
C...  No turbulence models, OR current field is not k or e.
C
      DO 1000 IEL=1,NEL
C
         IGRP=IGROUP(IEL)

         if (ifuservp) then
C
C           User specified fortran function   (pff 2/13/01)
            CALL NEKUVP (IEL)
            DIFMIN = VLMIN(VDIFF(1,1,1,IEL,IFIELD),NXYZ1)
            IF (DIFMIN .LE. 0.0) THEN
               WRITE (6,100) DIFMIN,IFIELD,IGRP
               CALL EXITT
            ENDIF
C
         ELSE IF(MATYPE(IGRP,IFIELD).EQ.1)THEN
C
C           Constant property within groups of elements
C
            CDIFF  = CPGRP(IGRP,IFIELD,1)
            CTRANS = CPGRP(IGRP,IFIELD,2)
            CALL CFILL(VDIFF (1,1,1,IEL,IFIELD),CDIFF,NXYZ1)
            CALL CFILL(VTRANS(1,1,1,IEL,IFIELD),CTRANS,NXYZ1)
            IF (CDIFF.LE.0.0) THEN
               WRITE(6,100) CDIFF,IFIELD,IGRP
  100          FORMAT(2X,'ERROR:  Non-positive diffusivity ('
     $        ,G12.3,') specified for field',I2,', group',I2
     $        ,' element',I4,'.'
     $        ,/,'ABORTING in VPROPS',//)
               CALL EXITT
            ENDIF
C
         ELSE IF(MATYPE(IGRP,IFIELD).EQ.2)THEN
C
C           User specified fortran function
C
            CALL NEKUVP (IEL)
C
            DIFMIN = VLMIN(VDIFF(1,1,1,IEL,IFIELD),NXYZ1)
            IF (DIFMIN .LE. 0.0) THEN
               WRITE (6,100) DIFMIN,IFIELD,IGRP
               CALL EXITT
            ENDIF
C
         ELSE IF(MATYPE(IGRP,IFIELD).EQ.0)THEN
C
C           Default constant property
C
            CDIFF  = CPFLD(IFIELD,1)
            CTRANS = CPFLD(IFIELD,2)
c           write(6,*) 'vdiff:',ifield,cdiff,ctrans
            CALL CFILL(VDIFF (1,1,1,IEL,IFIELD),CDIFF,NXYZ1)
            CALL CFILL(VTRANS(1,1,1,IEL,IFIELD),CTRANS,NXYZ1)
            IF (CDIFF.LE.0.0) THEN
               WRITE(6,200) CDIFF,IFIELD
  200          FORMAT(2X,'ERROR:  Non-positive diffusivity ('
     $        ,G12.3,') specified for field',I2,'.',/
     $        ,'ABORTING in VPROPS',//)
               CALL EXITT
            ENDIF
         ENDIF
C
 1000 CONTINUE
C
C     Turbulence models --- sum eddy viscosity/diffusivity
C
      IF (IFMODEL .AND. (IFIELD.EQ.1 .OR. IFIELD.EQ.2)) 
     $    CALL TVISCOS
C
      RETURN
      END
C
      SUBROUTINE NEKUVP (IEL)
C------------------------------------------------------------------
C
C     Generate user-specified material properties
C
C------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'PARALLEL'
      INCLUDE 'NEKUSE'
      ielg = lglel(iel)
c     IF (IFSTRS .AND. IFIELD.EQ.1) CALL STNRINV ! don't call! pff, 2007
      DO 10 K=1,NZ1
      DO 10 J=1,NY1
      DO 10 I=1,NX1
         CALL NEKASGN (I,J,K,IEL)
         CALL USERVP  (I,J,K,IELG)
         VDIFF (I,J,K,IEL,IFIELD) = UDIFF
         VTRANS(I,J,K,IEL,IFIELD) = UTRANS
 10   CONTINUE
      RETURN
      END
C
      SUBROUTINE DIAGNOS
      RETURN
      END
C
      SUBROUTINE SETSOLV
      INCLUDE 'SIZE'
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      IFSOLV = .FALSE.
      RETURN
      END

      SUBROUTINE MGGO
      RETURN
      END
      SUBROUTINE MGINIT
      RETURN
      END
C-----------------------------------------------------------------------
C
C     New files from LH on stress formulation (SFSUBS.FOR)
C
C-----------------------------------------------------------------------
      SUBROUTINE HMHZSF (NAME,U1,U2,U3,R1,R2,R3,H1,H2,
     $                   RMASK1,RMASK2,RMASK3,RMULT,
     $                   TOL,MAXIT,MATMOD)
C-----------------------------------------------------------------------
C
C     Compute solution to coupled Helmholtz equations 
C     (stress formulation)
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'TSTEP'
C
      DIMENSION U1(LX1,LY1,LZ1,1)
     $        , U2(LX1,LY1,LZ1,1)
     $        , U3(LX1,LY1,LZ1,1)
     $        , R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , R3(LX1,LY1,LZ1,1)
     $        , H1(LX1,LY1,LZ1,1)
     $        , H2(LX1,LY1,LZ1,1)
     $        , RMASK1(LX1,LY1,LZ1,1)
     $        , RMASK2(LX1,LY1,LZ1,1)
     $        , RMASK3(LX1,LY1,LZ1,1)
     $        , RMULT (LX1,LY1,LZ1,1)
      CHARACTER NAME*4
C
      NEL   = NELFLD(IFIELD)
      VOL   = VOLFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
c     sm1 = glsum(rmask1,ntot1)
c     sm2 = glsum(rmask2,ntot1)
c     sm3 = glsum(rmask3,ntot1)
c     write(6,*) ifield,ntot1,sm1,sm2,sm3,' SMASK'
c     write(6,*) ifield,ifldmhd
C
C     Current version uses only conjugate gradient iteration
C
C     IF ( .NOT.IFMGRID .OR. NAME.EQ.'NOMG') THEN
C

      CALL RMASK   (R1,R2,R3,NEL)
      CALL OPDSSUM (R1,R2,R3)

      CALL RZERO3  (U1,U2,U3,NTOT1)
C
      IF (IMESH.EQ.1) THEN
         CALL CHKTCGS (R1,R2,R3,RMASK1,RMASK2,RMASK3,RMULT,BINVM1,
     $                 VOL,TOL,NEL)
         CALL CGGOSF  (U1,U2,U3,R1,R2,R3,H1,H2,RMULT,BINVM1,
     $                 VOL,TOL,MAXIT,MATMOD)
      ELSE
         CALL CHKTCGS (R1,R2,R3,RMASK1,RMASK2,RMASK3,RMULT,BINTM1,
     $                 VOL,TOL,NEL)
         CALL CGGOSF  (U1,U2,U3,R1,R2,R3,H1,H2,RMULT,BINTM1,
     $                 VOL,TOL,MAXIT,MATMOD)
      ENDIF
C
C     ENDIF
C
      RETURN
      END
      SUBROUTINE CHKTCGS (R1,R2,R3,RMASK1,RMASK2,RMASK3,RMULT,BINV,
     $                    VOL,TOL,NEL)
C-------------------------------------------------------------------
C
C     Check that the tolerances are not too small for the CG-solver.
C     Important when calling the CG-solver (Gauss-Lobatto mesh) with
C     zero Neumann b.c.
C
C-------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'EIGEN'
      COMMON /CPRINT/ IFPRINT
      LOGICAL         IFPRINT
      COMMON /CTMP0/ WA (LX1,LY1,LZ1,LELT)
C
      DIMENSION R1    (LX1,LY1,LZ1,1)
     $        , R2    (LX1,LY1,LZ1,1)
     $        , R3    (LX1,LY1,LZ1,1)
     $        , RMASK1(LX1,LY1,LZ1,1)
     $        , RMASK2(LX1,LY1,LZ1,1)
     $        , RMASK3(LX1,LY1,LZ1,1)
     $        , RMULT (LX1,LY1,LZ1,1)
     $        , BINV  (LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
      IF (EIGAA .NE. 0.0) THEN
         ACONDNO = EIGGA/EIGAA
      ELSE
         ACONDNO = 10.0
      ENDIF
C
C     Check Single or double precision
C
      DELTA = 1.0E-9
      X     = 1.0 + DELTA
      Y     = 1.0
      DIFF  = ABS(X - Y)
      IF (DIFF .EQ. 0.0) EPS = 1.0E-6
      IF (DIFF .GT. 0.0) EPS = 1.0E-13
C
      CALL OPDOT (WA,R1,R2,R3,R1,R2,R3,NTOT1)
      RINIT = GLSC3(WA,BINV,RMULT,NTOT1)
      RINIT = SQRT (RINIT/VOL)
      RMIN  = EPS*RINIT
C
      IF (TOL.LT.RMIN) THEN
       TOLOLD = TOL
       TOL = RMIN
       IF (NID.EQ.0 .AND. IFPRINT)
     $ WRITE(6,*)'New CG1(stress)-tolerance (RINIT*epsm) = ',TOL,TOLOLD
      ENDIF
C
      IF (NDIM.EQ.2) THEN
         CALL ADD3 (WA,RMASK1,RMASK2,NTOT1)
      ELSE
         CALL ADD4 (WA,RMASK1,RMASK2,RMASK3,NTOT1)
      ENDIF
      BCNEU1 = GLSC2 (WA,RMULT,NTOT1)
      BCNEU2 = (NDIM) * GLSUM(RMULT,NTOT1)
      BCTEST = ABS(BCNEU1 - BCNEU2)
      IF (BCTEST .LT. 0.1) THEN
         IF (NDIM.EQ.2) THEN
            CALL ADD3 (WA,R1,R2,NTOT1)
         ELSE
            CALL ADD4 (WA,R1,R2,R3,NTOT1)
         ENDIF
         OTR    = GLSC2(WA,RMULT,NTOT1) / ( NDIM )
         TOLMIN = ABS(OTR) * ACONDNO
         IF (TOL .LT. TOLMIN) THEN
            TOLOLD = TOL
            TOL = TOLMIN
            IF (NID.EQ.0)
     $      WRITE (6,*) 'New CG1(stress)-tolerance (OTR) = ',TOL,TOLOLD
         ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE CGGOSF (U1,U2,U3,R1,R2,R3,H1,H2,RMULT,BINV,
     $                   VOL,tin,MAXIT,MATMOD)
C-----------------------------------------------------------------------
C
C     Conjugate gradient iteration for solution of coupled 
C     Helmholtz equations 
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      COMMON /SCREV/  DPC(LX1,LY1,LZ1,LELT)
     $     ,          P1 (LX1,LY1,LZ1,LELT)
      COMMON /SCRCH/  P2 (LX1,LY1,LZ1,LELT)
     $     ,          P3 (LX1,LY1,LZ1,LELT)
      COMMON /SCRMG/  PP1(LX1,LY1,LZ1,LELT)
     $     ,          PP2(LX1,LY1,LZ1,LELT)
     $     ,          PP3(LX1,LY1,LZ1,LELT)
     $     ,          WA (LX1,LY1,LZ1,LELT)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      COMMON /CPRINT/ IFPRINT
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV, IFPRINT
C
      DIMENSION U1   (LX1,LY1,LZ1,1)
     $        , U2   (LX1,LY1,LZ1,1)
     $        , U3   (LX1,LY1,LZ1,1)
     $        , R1   (LX1,LY1,LZ1,1)
     $        , R2   (LX1,LY1,LZ1,1)
     $        , R3   (LX1,LY1,LZ1,1)
     $        , H1   (LX1,LY1,LZ1,1)
     $        , H2   (LX1,LY1,LZ1,1)
     $        , RMULT(LX1,LY1,LZ1,1)
     $        , BINV (LX1,LY1,LZ1,1)
     $        , AP1  (LX1,LY1,LZ1,1)
     $        , AP2  (LX1,LY1,LZ1,1)
     $        , AP3  (LX1,LY1,LZ1,1)
      EQUIVALENCE (AP1,PP1),(AP2,PP2),(AP3,PP3)
C
      tol=tin
      if (param(22).ne.0) tol=abs(param(22))
      if (matmod.lt.0) tol = 1.e-4
c     if (matmod.lt.0) tol = tin
c
      NEL   = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
C
C     Set logical flags
C
      IF ( .NOT.IFSOLV ) THEN
           CALL SETFAST (H1,H2,IMESH)
           IFSOLV = .TRUE.
      ENDIF
      CALL OPDOT (WA,R1,R2,R3,R1,R2,R3,NTOT1)
      RBNORM = GLSC3(WA,BINV,RMULT,NTOT1)
      RBNORM = SQRT ( RBNORM / VOL )
      IF (RBNORM .LT. TOL) THEN
C         IF ( .NOT.IFPRINT )  GOTO 9999
          IF (MATMOD.GE.0.AND.NID.EQ.0) WRITE (6,2000) RBNORM
          IF (MATMOD.LT.0.AND.NID.EQ.0) WRITE (6,2010) RBNORM
          GOTO 9999
      ENDIF
C
C     Evaluate diagional pre-conidtioner for fluid solve
C
c     write(6,*) matmod,ifield,ifldmhd,'  solver'
      IF (MATMOD.GE.0) THEN
          CALL SETPREC (DPC,H1,H2,IMESH,1)
          CALL SETPREC (WA ,H1,H2,IMESH,2)
          CALL ADD2    (DPC,WA,NTOT1)
          IF (NDIM.EQ.3) THEN
             CALL SETPREC (WA,H1,H2,IMESH,3)
             CALL ADD2    (DPC,WA,NTOT1)
          ENDIF
      ELSE
          CALL RONE (DPC,NTOT1)
      ENDIF
C
      CALL COL3 (PP1,DPC,R1,NTOT1)
      CALL COL3 (PP2,DPC,R2,NTOT1)
      CALL COPY (P1,PP1,NTOT1)
      CALL COPY (P2,PP2,NTOT1)
      IF (NDIM.EQ.3) THEN
         CALL COL3 (PP3,DPC,R3,NTOT1)
         CALL COPY (P3,PP3,NTOT1)
      ENDIF
C
      CALL OPDOT  (WA,R1,R2,R3,PP1,PP2,PP3,NTOT1)
      RPP1 = GLSC2(WA,RMULT,NTOT1)
C
      DO 1000 ITER=1,MAXIT
C
         CALL AXHMSF  (AP1,AP2,AP3,P1,P2,P3,H1,H2,MATMOD)
         CALL RMASK   (AP1,AP2,AP3,NEL)
         CALL OPDSSUM (AP1,AP2,AP3)
C
         CALL OPDOT  (WA,P1,P2,P3,AP1,AP2,AP3,NTOT1)
         PAP   = GLSC2(WA,RMULT,NTOT1)
         ALPHA = RPP1 / PAP
         CALL OPADDS (U1,U2,U3,P1 ,P2 ,P3 , ALPHA,NTOT1,2)
         CALL OPADDS (R1,R2,R3,AP1,AP2,AP3,-ALPHA,NTOT1,2)
C
         CALL OPDOT  (WA,R1,R2,R3,R1,R2,R3,NTOT1)
         RBNORM = GLSC3(WA,BINV,RMULT,NTOT1)
         RBNORM = SQRT (RBNORM/VOL)
         if (iter.eq.1) r0 = rbnorm
c        if (matmod.lt.0) write(6,9) iter,rbnorm,r0
c  9     format(i5,1p2e14.5,' mesh')
         IF (RBNORM .LT. TOL) THEN
            IFIN = ITER - 1
C           IF ( .NOT.IFPRINT )  GOTO 9999
            if (nid.eq.0) then
               if (matmod.ge.0) write(6,3000) ifin,rbnorm,tol,r0
               if (matmod.lt.0) write(6,3010) ifin,rbnorm,tol,r0
            endif
            goto 9999
         ENDIF
C
         CALL COL3 (PP1,DPC,R1,NTOT1)
         CALL COL3 (PP2,DPC,R2,NTOT1)
         IF (NDIM.EQ.3) CALL COL3 (PP3,DPC,R3,NTOT1)
C
         CALL OPDOT (WA,R1,R2,R3,PP1,PP2,PP3,NTOT1)
         RPP2 = RPP1
         RPP1 = GLSC2(WA,RMULT,NTOT1)
         BETA = RPP1/RPP2
         CALL OPADDS (P1,P2,P3,PP1,PP2,PP3,BETA,NTOT1,1)
C
 1000 CONTINUE
C 
      IF (MATMOD.GE.0.AND.NID.EQ.0) WRITE (6,3001) ITER,RBNORM
      IF (MATMOD.LT.0.AND.NID.EQ.0) WRITE (6,3011) ITER,RBNORM
C
 9999 CONTINUE
      IFSOLV = .FALSE.
C
c     call outpost(u1,u2,u3,u3,u3,'   ')
c     call outpost(h1,h2,u3,u3,u3,'   ')
c     call exitt
c
      return
C
 2000 FORMAT(13X,'Helmholtz3/fluid: no iteration - RBNORM =', E13.4)
 2010 FORMAT(13X,'Helmholtz3/ mesh: no iteration - RBNORM =', E13.4)
 3000 FORMAT(13X,'Helmholtz3/fluid: ',I6,5E13.4)
 3010 FORMAT(13X,'Helmholtz3/ Mesh: ',I6,5E13.4)
 3001 FORMAT(I6,' Failed to converge in HMHZSF/Fluid: RBNORM =',E13.6)
 3011 FORMAT(I6,' Failed to converge in HMHZSF/Mesh : RBNORM =',E13.6)
C
      end

      SUBROUTINE AXHMSF (AU1,AU2,AU3,U1,U2,U3,H1,H2,MATMOD)
C-----------------------------------------------------------------------
C
C     Compute the coupled Helmholtz matrix-vector products
C
C     Fluid (MATMOD .GE. 0) :  Hij Uj = Aij*Uj + H2*B*Ui 
C     Solid (MATMOD .LT. 0) :  Hij Uj = Kij*Uj
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'TSTEP'
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      DIMENSION AU1(LX1,LY1,LZ1,1)
     $        , AU2(LX1,LY1,LZ1,1)
     $        , AU3(LX1,LY1,LZ1,1)
     $        , U1 (LX1,LY1,LZ1,1)
     $        , U2 (LX1,LY1,LZ1,1)
     $        , U3 (LX1,LY1,LZ1,1)
     $        , H1 (LX1,LY1,LZ1,1)
     $        , H2 (LX1,LY1,LZ1,1)
C

      if (ifsplit) then
         call axhelm(au1,u1,h1,h2,1,1)
         call axhelm(au2,u2,h1,h2,1,2)
         if (if3d) call axhelm(au3,u3,h1,h2,1,3)
         return
      endif

      NEL   = NELFLD(IFIELD)
      NTOT1 = NX1*NY1*NZ1*NEL
C
C     Calculate coupled  Aij Uj  products
C
      IF ( .NOT.IFSOLV ) CALL SETFAST (H1,H2,IMESH)
      CALL STNRATE (U1,U2,U3,NEL,MATMOD)
      CALL STRESS  (H1,H2,NEL,MATMOD,IFAXIS)
      CALL AIJUJ   (AU1,AU2,AU3,NEL,IFAXIS)
C
C     Add other Helmholtz contributions
C
      IF (IFH2 .AND. MATMOD.GE.0) THEN
         CALL ADDCOL4 (AU1,BM1,H2,U1,NTOT1)
         CALL ADDCOL4 (AU2,BM1,H2,U2,NTOT1)
         IF (NDIM.EQ.3) CALL ADDCOL4 (AU3,BM1,H2,U3,NTOT1)
      ENDIF
C
      RETURN
      END
      SUBROUTINE STNRATE (U1,U2,U3,NEL,MATMOD)
C
C     Compute strainrates
C
C     CAUTION : Stresses and strainrates share the same scratch commons
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      COMMON /CTMP0/ EXZ(LX1,LY1,LZ1,LELT)
     $             , EYZ(LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ EXX(LX1,LY1,LZ1,LELT)
     $             , EXY(LX1,LY1,LZ1,LELT)
     $             , EYY(LX1,LY1,LZ1,LELT)
     $             , EZZ(LX1,LY1,LZ1,LELT)
C
      DIMENSION U1(LX1,LY1,LZ1,1)
     $        , U2(LX1,LY1,LZ1,1)
     $        , U3(LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
      CALL RZERO3 (EXX,EYY,EZZ,NTOT1)
      CALL RZERO3 (EXY,EXZ,EYZ,NTOT1)
C
      CALL UXYZ  (U1,EXX,EXY,EXZ,NEL)
      CALL UXYZ  (U2,EXY,EYY,EYZ,NEL)
      IF (NDIM.EQ.3) CALL UXYZ   (U3,EXZ,EYZ,EZZ,NEL)
C
      CALL INVCOL2 (EXX,JACM1,NTOT1)
      CALL INVCOL2 (EXY,JACM1,NTOT1)
      CALL INVCOL2 (EYY,JACM1,NTOT1)
C
      IF (IFAXIS) CALL AXIEZZ (U2,EYY,EZZ,NEL)
C
      IF (NDIM.EQ.3) THEN
         CALL INVCOL2 (EXZ,JACM1,NTOT1)
         CALL INVCOL2 (EYZ,JACM1,NTOT1)
         CALL INVCOL2 (EZZ,JACM1,NTOT1)
      ENDIF
C
      RETURN
      END
      SUBROUTINE UXYZ (U,EX,EY,EZ,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      COMMON /SCRSF/ UR(LX1,LY1,LZ1,LELT)
     $             , US(LX1,LY1,LZ1,LELT)
     $             , UT(LX1,LY1,LZ1,LELT)
C
      DIMENSION U (LX1,LY1,LZ1,1)
     $        , EX(LX1,LY1,LZ1,1)
     $        , EY(LX1,LY1,LZ1,1)
     $        , EZ(LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
      CALL URST (U,UR,US,UT,NEL)
C
      CALL ADDCOL3 (EX,RXM1,UR,NTOT1)
      CALL ADDCOL3 (EX,SXM1,US,NTOT1)
      CALL ADDCOL3 (EY,RYM1,UR,NTOT1)
      CALL ADDCOL3 (EY,SYM1,US,NTOT1)
C
      IF (NDIM.EQ.3) THEN
         CALL ADDCOL3 (EZ,RZM1,UR,NTOT1)
         CALL ADDCOL3 (EZ,SZM1,US,NTOT1)
         CALL ADDCOL3 (EZ,TZM1,UT,NTOT1)
         CALL ADDCOL3 (EX,TXM1,UT,NTOT1)
         CALL ADDCOL3 (EY,TYM1,UT,NTOT1)
      ENDIF
C
      RETURN
      END
      SUBROUTINE URST (U,UR,US,UT,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
C
      DIMENSION U (LX1,LY1,LZ1,1)
     $        , UR(LX1,LY1,LZ1,1)
     $        , US(LX1,LY1,LZ1,1)
     $        , UT(LX1,LY1,LZ1,1)
C
      DO 100 IEL=1,NEL
         IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
         CALL DDRST (U (1,1,1,IEL),UR(1,1,1,IEL),
     $               US(1,1,1,IEL),UT(1,1,1,IEL))
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE DDRST (U,UR,US,UT)
C
      INCLUDE 'SIZE'
      INCLUDE 'DXYZ'
C
      DIMENSION U (LX1,LY1,LZ1)
     $        , UR(LX1,LY1,LZ1)
     $        , US(LX1,LY1,LZ1)
     $        , UT(LX1,LY1,LZ1)
C
      NXY1 = NX1*NY1
      NYZ1 = NY1*NZ1
C
      CALL MXM (DXM1,NX1,U,NX1,UR,NYZ1)
      IF (NDIM.EQ.2) THEN
         CALL MXM (U,NX1,DYTM1,NY1,US,NY1)
      ELSE
         DO 10 IZ=1,NZ1
         CALL MXM (U(1,1,IZ),NX1,DYTM1,NY1,US(1,1,IZ),NY1)
   10    CONTINUE
         CALL MXM (U,NXY1,DZTM1,NZ1,UT,NZ1)
      ENDIF
C
      RETURN
      END
      SUBROUTINE AXIEZZ (U2,EYY,EZZ,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
C
      DIMENSION U2 (LX1,LY1,LZ1,1)
     $        , EYY(LX1,LY1,LZ1,1)
     $        , EZZ(LX1,LY1,LZ1,1)
C
      NXYZ1  = NX1*NY1*NZ1
C
      DO 100 IEL=1,NEL
         IF ( IFRZER(IEL) ) THEN
            DO 200 IX=1,NX1
               EZZ(IX, 1,1,IEL) = EYY(IX,1,1,IEL)
            DO 200 IY=2,NY1
               EZZ(IX,IY,1,IEL) = U2(IX,IY,1,IEL) / YM1(IX,IY,1,IEL)
  200       CONTINUE
         ELSE
            CALL INVCOL3 (EZZ(1,1,1,IEL),U2(1,1,1,IEL),YM1(1,1,1,IEL),
     $                    NXYZ1)
         ENDIF
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine flush_io
      return
      end
c-----------------------------------------------------------------------
      subroutine fcsum2(xsum,asum,x,e,f)
c
c     Compute the weighted sum of X over face f of element e
c
c     x is an (NX,NY,NZ) data structure
c     f  is in the preprocessor notation 
c
c     xsum is sum (X*area)
c     asum is sum (area)


      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'
      real x(lx1,ly1,lz1,1)
      integer e,f,fd

      asum = 0.
      xsum = 0.

c     Set up counters ;  fd is the dssum notation.
      call dsset(nx1,ny1,nz1)
      fd     = eface1(f)
      js1    = skpdat(1,fd)
      jf1    = skpdat(2,fd)
      jskip1 = skpdat(3,fd)
      js2    = skpdat(4,fd)
      jf2    = skpdat(5,fd)
      jskip2 = skpdat(6,fd)

      i = 0
      do j2=js2,jf2,jskip2
      do j1=js1,jf1,jskip1
         i = i+1
         xsum = xsum+area(i,1,f,e)*x(j1,j2,1,e)
         asum = asum+area(i,1,f,e)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      function surf_mean(u,ifld,bc_in,ierr)

      include 'SIZE'
      include 'TOTAL'

      real u(1)

      integer e,f
      character*3 bc_in

      usum = 0
      asum = 0

      nface = 2*ndim
      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,ifld).eq.bc_in) then
            call fcsum2(usum_f,asum_f,u,e,f)
            usum = usum + usum_f
            asum = asum + asum_f
         endif
      enddo
      enddo

      usum = glsum(usum,1)  ! sum across processors
      asum = glsum(asum,1)

      surf_mean = usum
      ierr      = 1

      if (asum.gt.0) surf_mean = usum/asum
      if (asum.gt.0) ierr      = 0

      return
      end
c-----------------------------------------------------------------------
