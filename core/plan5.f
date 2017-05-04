c-----------------------------------------------------------------------
      subroutine plan5(igeom)

c     Two-step Richardson Extrapolation.
c     Operator splitting technique.

      include 'SIZE'
      include 'TOTAL'

      common /scrns/  resv  (lx1*ly1*lz1*lelv,3)

      n   = nx1*ny1*nz1*nelv
      n2  = nx2*ny2*nz2*nelv
      dt2 = dt/2
      dti = 1/dt

      if (igeom.eq.2) then

      if (ifmvbd) call opcopy
     $  (wxlag(1,1,1,1,2),wylag(1,1,1,1,2),wzlag(1,1,1,1,2),xm1,ym1,zm1)
      do i=1,n
         s = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti  ! Add  density*mass/dt,
         vxlag(i,1,1,1,2)=s*vx(i,1,1,1)          ! equivalent to using
         vylag(i,1,1,1,2)=s*vy(i,1,1,1)          ! density*mass/(dt/2)
         vzlag(i,1,1,1,2)=s*vz(i,1,1,1)          ! in the first place...
      enddo

      call pn2_step(vxlag,vylag,vzlag,prlag,0,dt)  ! One step of Pn-Pn-2

      do i=1,n                                          ! Add  density*mass/dt,
         bfx(i,1,1,1)=bfx(i,1,1,1)+vxlag(i,1,1,1,2)     ! equivalent to using
         bfy(i,1,1,1)=bfy(i,1,1,1)+vylag(i,1,1,1,2)     ! density*mass/(dt/2)
         bfz(i,1,1,1)=bfz(i,1,1,1)+vzlag(i,1,1,1,2)     ! in the first place...
      enddo

      if (ifmvbd) then
         write (*,*) 'ifmbvd is true'
        call opcopy
     $  (xm1,ym1,zm1,wxlag(1,1,1,1,2),wylag(1,1,1,1,2),wzlag(1,1,1,1,2))
        call geom_reset(0)
      endif

      time = time-dt2
      call pn2_step(vx,vy,vz,pr,1,dt2)      ! One step of Pn-Pn-2, dt/2

      time = time+dt2
      call setup_convect(2)  ! Map vx --> vxd
      call setprop

      call pn2_step(vx,vy,vz,pr,0,dt2)      ! One step of Pn-Pn-2, dt/2

      do i=1,n
         vx(i,1,1,1)=2*vx(i,1,1,1)-vxlag(i,1,1,1,1)
         vy(i,1,1,1)=2*vy(i,1,1,1)-vylag(i,1,1,1,1)
         vz(i,1,1,1)=2*vz(i,1,1,1)-vzlag(i,1,1,1,1)
      enddo

      do i=1,n2
         pr(i,1,1,1)=2*pr(i,1,1,1)-prlag(i,1,1,1,1)
      enddo

      call ortho(pr)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pn2_step(ux,uy,uz,pu,iresv,dtl)
      include 'SIZE'
      include 'TOTAL'

      parameter (lv=lx1*ly1*lz1*lelt)
      real ux(1),uy(1),uz(1),pu(1)

      common /scrns/  resv  (lx1*ly1*lz1*lelv,3)
     $ ,              dv1   (lx1*ly1*lz1*lelv)
     $ ,              dv2   (lx1*ly1*lz1*lelv)
     $ ,              dv3   (lx1*ly1*lz1*lelv)
      common /scrvh/  h1    (lx1*ly1*lz1*lelv)
     $ ,              h2    (lx1*ly1*lz1*lelv)


      if (nx1.eq.nx2) 
     $   call exitti('pn2_step requires lx2=lx1-2 in SIZE$',nx2)

      ifield = 1                ! Set field for velocity
      n   = nx1*ny1*nz1*nelv
      n2  = nx2*ny2*nz2*nelv

      dti = 1/dtl
      call copy    (h1,vdiff ,n)
      call cmult2  (h2,vtrans,dti,n)

      if (iresv.eq.0) then ! bfx etc is preserved if iresv=1

                    call makeuf  ! Initialize bfx, bfy, bfz
        if (ifmvbd) call admeshv ! Add div(W.u_i)

        call convop(resv(1,1),vx)
        call convop(resv(1,2),vy)
        call convop(resv(1,3),vz)

        do i=1,n
           b=vtrans(i,1,1,1,1)*bm1(i,1,1,1)
           s=1./dtl
           bfx(i,1,1,1)=bfx(i,1,1,1)+b*(s*vx(i,1,1,1)-resv(i,1))
           bfy(i,1,1,1)=bfy(i,1,1,1)+b*(s*vy(i,1,1,1)-resv(i,2))
           bfz(i,1,1,1)=bfz(i,1,1,1)+b*(s*vz(i,1,1,1)-resv(i,3))
        enddo

      endif

      call adv_geom(dtl) ! Advance the geometry

      call opcopy  (ux,uy,uz,vx,vy,vz)
      call bcdirvc (ux,uy,uz,v1mask,v2mask,v3mask) ! Don't forget bcneutr !
      call ophx    (resv(1,1),resv(1,2),resv(1,3),ux,uy,uz,h1,h2)

      do i=1,n
         resv(i,1)=bfx(i,1,1,1)-resv(i,1) ! rhs = rhs - H*u
         resv(i,2)=bfy(i,1,1,1)-resv(i,2)
         resv(i,3)=bfz(i,1,1,1)-resv(i,3)
      enddo

      tolhv = abs(param(22))
      nmxh  = 100
      call ophinv_pr(dv1,dv2,dv3
     $   ,resv(1,1),resv(1,2),resv(1,3),h1,h2,tolhv,nmxh)

      call opadd2  (ux,uy,uz,dv1,dv2,dv3)

      bd(1) = 1.0
      call rzero(pu,n2)

      dt_old = dt
      dt = dtl
      call incomprn(ux,uy,uz,pu)
      dt = dt_old

      return
      end
c-----------------------------------------------------------------------
      subroutine adv_geom(dtl) ! Advance the geometry
      include 'SIZE'
      include 'TOTAL'

      param(28) = 1   !  This forces Euler Forward for Mesh Update
                      !  Note: p28 must be set prior to call settime

      dt_tmp = dt     ! Save "full" dt value
      dt     = dtl

      ifield = 1
      call gengeomo1(2)
c     call gengeom(2)
      ifield = 1

      dt     = dt_tmp ! Restore dt

      return
      end
c-----------------------------------------------------------------------
      subroutine gengeomo1(igeom)
C----------------------------------------------------------------------
C
C     Generate geometry data
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
      include 'WZ'
C
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C

      if (nio.eq.0.and.istep.le.1) write(6,*) 'generate geometry data'

      IF (IGEOM.EQ.1) THEN
         RETURN
      ELSEIF (IGEOM.EQ.2) THEN
         CALL LAGMASS
         IF (ISTEP.EQ.0) CALL GENCOOR (XM3,YM3,ZM3)
         if (istep.ge.1) call updcoor1
         CALL GEOM1 (XM3,YM3,ZM3)
         CALL GEOM2
         CALL UPDMSYS (1)
         CALL VOLUME
         CALL SETINVM
         CALL SETDEF
         CALL SFASTAX
         IF (ISTEP.GE.1) CALL EINIT
      ELSEIF (IGEOM.EQ.3) THEN
c
c        Take direct stiffness avg of mesh
c
         ifieldo = ifield
         if (.not.ifneknekm) CALL GENCOOR (XM3,YM3,ZM3)
         if (ifheat) then
            ifield = 2
            CALL dssum(xm3,nx3,ny3,nz3)
            call col2 (xm3,tmult,ntot3)
            CALL dssum(ym3,nx3,ny3,nz3)
            call col2 (ym3,tmult,ntot3)
            if (if3d) then
               CALL dssum(xm3,nx3,ny3,nz3)
               call col2 (xm3,tmult,ntot3)
            endif
         else
            ifield = 1
            CALL dssum(xm3,nx3,ny3,nz3)
            call col2 (xm3,vmult,ntot3)
            CALL dssum(ym3,nx3,ny3,nz3)
            call col2 (ym3,vmult,ntot3)
            if (if3d) then
               CALL dssum(xm3,nx3,ny3,nz3)
               call col2 (xm3,vmult,ntot3)
            endif
         endif
         CALL GEOM1 (XM3,YM3,ZM3)
         CALL GEOM2
         CALL UPDMSYS (1)
         CALL VOLUME
         CALL SETINVM
         CALL SETDEF
         CALL SFASTAX
         ifield = ifieldo
      ENDIF

      if (nio.eq.0.and.istep.le.1) then
        write(6,*) 'done :: generate geometry data' 
        write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine updcoor1
C-----------------------------------------------------------------------
C
C     Subroutine to update geometry for moving boundary problems
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TSTEP'
C
      IFIELD = 0
      NEL    = NELFLD(IFIELD)
C
C     Update collocation points coordinates
C
      call updxyz1(nel)
C
C     Shift lagged mesh velocity
C
c     CALL LAGMSHV (NEL)
C
      return
      end
C-----------------------------------------------------------------------
      subroutine updxyz1(nel)
C
      include 'SIZE'
      include 'TSTEP'
      include 'MVGEOM'
      include 'GEOM'
      COMMON /SCRSF/ UX(LX1,LY1,LZ1,LELT)
     $             , UY(LX1,LY1,LZ1,LELT)
     $             , UZ(LX1,LY1,LZ1,LELT)
      DIMENSION ABM(3)
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
      write (*,*) 'dt,updxyz=',dt
      DO 10 I=1,NBD
   10 ABM(I) = DT*ABMSH(I)
C
      IF (ISTEP.EQ.0) THEN
         CALL COPY (UX,WX,NTOT1)
         CALL COPY (UY,WY,NTOT1)
         IF (NDIM.EQ.3) CALL COPY (UZ,WZ,NTOT1)
      ELSE
         call cmult2 (ux,wx,dt,ntot1)
         call cmult2 (uy,wy,dt,ntot1)
         if (ndim.eq.3) call cmult2 (uz,wz,dt,ntot1)
c        DO 100 ILAG=2,NBD
c           CALL ADD2S2 (UX,WXLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
c           CALL ADD2S2 (UY,WYLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
c           IF (NDIM.EQ.3) 
c    $      CALL ADD2S2 (UZ,WZLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
c 100    CONTINUE
      ENDIF
C
      CALL ADD2 (XM1,UX,NTOT1)
      CALL ADD2 (YM1,UY,NTOT1)
      IF (NDIM.EQ.3) CALL ADD2 (ZM1,UZ,NTOT1)
C
      return
      end
C-----------------------------------------------------------------------
