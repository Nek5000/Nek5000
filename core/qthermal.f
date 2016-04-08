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

      common /scrns/ w1(lx1,ly1,lz1,lelt)
     $              ,w2(lx1,ly1,lz1,lelt)
     $              ,tx(lx1,ly1,lz1,lelt)
     $              ,ty(lx1,ly1,lz1,lelt)
     $              ,tz(lx1,ly1,lz1,lelt)

      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelv

      if (.not.iflomach) then
         call rzero(qtl,ntot)
         return
      endif

      ifld_save = ifield

c - - Assemble RHS of T-eqn
      ifield=2
      call setqvol (qtl) ! volumetric heating source
      call col2    (qtl,bm1,ntot)

      ifield=1     !set right gs handle (QTL is only defined on the velocity mesh)
      call opgrad  (tx,ty,tz,t)
      call opdssum (tx,ty,tz)
      call opcolv  (tx,ty,tz,binvm1)
      call opcolv  (tx,ty,tz,vdiff(1,1,1,1,2))
      call opdiv   (w2,tx,ty,tz)

      call add2    (qtl,w2,ntot)
      call dssum   (qtl,nx1,ny1,nz1)
      call col2    (qtl,binvm1,ntot)

      ! QTL = T_RHS/(rho*cp**T)

      call col3    (w2,vtrans(1,1,1,1,2),t,ntot)
      call invcol2 (qtl,w2,ntot)

      if (ifcvode .and. ifvcor) then
         call make_p0th
         call add_qthermal
      else
         dp0thdt = 0.0
      endif

      ifield = ifld_save

      return
      end
