c-------------------------------------------------------------------------
      subroutine qthermal

C     Compute the thermal divergence QTL 
C
C     QTL := div(v) = -1/rho * Drho/Dt
c
c     If we use the ideal gas law and assume
c     that p,R is const we end up with
c     QTL = 1/(rho*cp) rho*cp*DT/Dt
C
C     where rho*cp*DT/Dt represents the RHS of the
C     energy equation expressed in terms of temperature.

      include 'SIZE'
      include 'TOTAL'

      COMMON /SCRNS/ w1(LX1,LY1,LZ1,LELT)
     $              ,w2(LX1,LY1,LZ1,LELT)
     $              ,tx(LX1,LY1,LZ1,LELT)
     $              ,ty(LX1,LY1,LZ1,LELT)
     $              ,tz(LX1,LY1,LZ1,LELT)

      ntot = nx1*ny1*nz1*nelv

      if (.not.iflomach) then
         call rzero(qtl,ntot)
         return
      endif

      ifld_save = ifield

c - - Assemble RHS of T-eqn
      ifield=2
      call setqvol (QTL) ! volumetric heating source
      call col2    (QTL,BM1,ntot)

      ifield=1     !set right gs handle (QTL is only defined on the velocity mesh)
      call opgrad  (tx,ty,tz,T)
      call opdssum (tx,ty,tz)
      call opcolv  (tx,ty,tz,binvm1)
      call opcolv  (tx,ty,tz,vdiff(1,1,1,1,2))
      call opdiv   (w2,tx,ty,tz)

      call add2    (QTL,w2,ntot)

      call col3    (w2,vtrans(1,1,1,1,2),T,ntot)
      call invcol2 (QTL,w2,ntot)

      call dssum   (QTL,nx1,ny1,nz1)
      call col2    (QTL,binvm1,ntot)

      ifield = ifld_save

      return
      end
