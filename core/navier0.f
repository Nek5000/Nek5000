      SUBROUTINE ESOLVER (RES,H1,H2,H2INV,INTYPE)
C---------------------------------------------------------------------
C
C     Choose E-solver
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'ESOLV'
      INCLUDE 'INPUT'
C
      REAL RES   (LX2,LY2,LZ2,LELV)
      REAL H1    (LX1,LY1,LZ1,LELV)
      REAL H2    (LX1,LY1,LZ1,LELV)
      REAL H2INV (LX1,LY1,LZ1,LELV)
      common /scruz/ wk1(lx2*ly2*lz2*lelv)
     $             , wk2(lx2*ly2*lz2*lelv)
     $             , wk3(lx2*ly2*lz2*lelv)
c
      include 'CTIMER'
c
      if (icalld.eq.0) teslv=0.0
      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()
C
c     write(6,*) solver_type,' solver type',iesolv
      IF (IESOLV.EQ.1) THEN
         if (solver_type.eq.'fdm') then
            ntot2 = nx2*ny2*nz2*nelv
            call gfdm_pres_solv  (wk1,res,wk2,wk3)
            call copy            (res,wk1,ntot2)
         else
            if (param(42).eq.1.or.solver_type.eq.'pdm') then
               CALL UZAWA (RES,H1,H2,H2INV,INTYPE,ICG)
            else
               call uzawa_gmres(res,h1,h2,h2inv,intype,icg)
            endif
         endif
      ELSE
         WRITE(6,*) 'ERROR: E-solver does not exist',iesolv
         WRITE(6,*) 'Stop in ESOLVER'
         CALL EXITT
      ENDIF
      teslv=teslv+(dnekclock()-etime1)
      RETURN
      END
      SUBROUTINE ESTRAT
C---------------------------------------------------------------------------
C
C     Decide strategy for E-solver
C
C---------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
C
      IESOLV = 1
      if (ifsplit) iesolv=0
      IF (ISTEP.LT.2.AND.NID.EQ.0) WRITE(6,10) IESOLV
   10 FORMAT(2X,'E-solver strategy:',I2)

c
      solver_type='itr'
      if (param(116).ne.0) solver_type='fdm'
c     if (param(90).ne.0)  solver_type='itn'
c
c     The following change recognizes that geometry is logically 
c     tensor-product, but deformed:  pdm = Preconditioner is fdm
c
      if (param(59).ne.0.and.solver_type.eq.'fdm') solver_type='pdm'
C
      RETURN
      END
      SUBROUTINE EINIT
C-----------------------------------------------------------------------------
C
C     Initialize E-solver
C
C-----------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ESOLV'
      COMMON /SCRHI/  H2INV (LX1,LY1,LZ1,LELV)
      LOGICAL IFNEWP
C
      CALL ESTRAT 
      RETURN
      END
