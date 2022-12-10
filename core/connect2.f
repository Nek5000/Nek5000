c-----------------------------------------------------------------------
      subroutine readat

c     Read data from preprocessor input files (rea,par,re2,co2,ma2,etc.)

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'RESTART'

      logical parfound

      call flush_io

      ! parfle is in 'INPUT'
      if (nid.eq.0) inquire(file=parfle, exist=parfound)
      call bcast(parfound,lsize) ! check for par file

      get_vert_called = 0

      if (parfound) then
        if(nio.eq.0) write(6,'(a,a)') ' Reading ', parfle

        call setDefaultParam

        if (nid.eq.0) call par_read(ierr)
        call bcast(ierr, isize)
        if (ierr.ne.0) call exitt
        call bcastParam

        call usrdat0

        if (ifnewre2reader) then
          if(nio.eq.0) write(6,'(a)') ' Using new re2 reader ...'
          call readat_big_v2
        else
          if(nio.eq.0) write(6,'(a)') ' Using old re2 reader ...'
          call readat_par
        endif

      else
        if(nio.eq.0) write(6,'(a,a)') ' Reading .rea file '
        call readat_std
      endif

      call set_boundary_ids

      return
      end
c-----------------------------------------------------------------------
      subroutine readat_std
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'RESTART'

      logical ifbswap,ifre2
      character*132 string
      integer idum(3*numsts+3)

      etime0 = dnekclock_sync()

      ierr = 0
      if(nid.eq.0) then
        write(6,'(A,A)') ' Reading ', reafle
        open (unit=9,file=reafle,status='old', iostat=ierr)
      endif

      call bcast(ierr,isize)
      if (ierr .gt. 0) call exitti('Cannot open rea file!$',1)

C     Read parameters and logical flags
      call rdparam

C     Read Mesh Info 
      if(nid.eq.0) then
        read(9,*)    ! xfac,yfac,xzero,yzero
        read(9,*)    ! dummy
        read(9,*)  nelgs,ldimr,nelgv
        nelgt = abs(nelgs)
      endif
      call bcast(ldimr,ISIZE)
      call bcast(nelgs,ISIZE)
      call bcast(nelgv,ISIZE)
      call bcast(nelgt,ISIZE)
      ifre2 = .false.
      if (nelgs.lt.0) ifre2 = .true.

      call usrdat0

      if (nelgt.gt.350000 .and. .not.ifre2) 
     $   call exitti('Problem size requires .re2!$',1)

      if (ifre2) call read_re2_hdr(ifbswap, .true.) ! rank0 will open and read
      call chk_nel  ! make certain sufficient array sizes

      call mapelpr

      if (ifre2) then
        call read_re2_data(ifbswap, .true., .true., .true.)
      else
        maxrd = 32               ! max # procs to read at once
        mread = (np-1)/maxrd+1   ! mod param
        iread = 0                ! mod param
        x     = 0
        do i=0,np-1,maxrd
           call nekgsync()
           if (mod(nid,mread).eq.iread) then
              if (nid.ne.0) then
                open(UNIT=9,FILE=REAFLE,STATUS='OLD')
                call cscan(string,'MESH DATA',9)
                read(9,*) string
              endif 
              call rdmesh
              call rdcurve !  Curved side data
              call rdbdry  !  Boundary Conditions
              if (nid.ne.0) close(unit=9)
           endif
           iread = iread + 1
        enddo
      endif

C     Read Restart options / Initial Conditions / Drive Force
      CALL RDICDF
C     Read materials property data
      CALL RDMATP
C     Read history data
      CALL RDHIST
C     Read output specs
      CALL RDOUT
C     Read objects
      CALL RDOBJ

      call nekgsync()

C     End of input data, close read file.
      if(nid.eq.0) then
        close(unit=9)
        call echopar
        write(6,'(A,g13.5,A,/)')  ' done :: read .rea file ',
     $                             dnekclock()-etime0,' sec'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_boundary_ids

      include 'SIZE'
      include 'TOTAL'

      call izero(boundaryID,  size(boundaryID))
      call izero(boundaryIDt, size(boundaryIDt))

      ifld = 2 
      if(ifflow) ifld = 1
      do iel = 1,nelv
      do ifc = 1,2*ldim   
         boundaryID(ifc,iel) = bc(5,ifc,iel,ifld)
      enddo
      enddo

      ntmsh = 0
      do i=1,ldimt
         if(iftmsh(1+i)) ntmsh = ntmsh + 1 
      enddo

      if (ntmsh.gt.0) then
        do iel = 1,nelt
        do ifc = 1,2*ldim   
           boundaryIDt(ifc,iel) = bc(5,ifc,iel,2)
        enddo
        enddo
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine readat_old
C
C     Read in data from preprocessor input file (.rea)
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'
      include 'CTIMER'
 
      logical ifbswap,ifre2,parfound
      character*132 string
      integer idum(3*numsts+3)

      ierr = 0
      call flush_io

      ! check if new rea file version exists
      if(nid.eq.0) inquire(file=parfle, exist=parfound)
      call bcast(parfound,lsize)
      if (parfound) then
         if(nio.eq.0) write(6,'(A,A)') ' Reading ', parfle
         call readat_par
         goto 99
      endif  

      etime0 = dnekclock_sync()

      if(nid.eq.0) then
        write(6,'(A,A)') ' Reading ', reafle
        open (unit=9,file=reafle,status='old', iostat=ierr)
      endif

      call bcast(ierr,isize)
      if (ierr .gt. 0) call exitti('Cannot open rea file!$',1)

C     Read parameters and logical flags
      call rdparam

C     Read Mesh Info 
      if(nid.eq.0) then
        read(9,*)    ! xfac,yfac,xzero,yzero
        read(9,*)    ! dummy
        read(9,*)  nelgs,ldimr,nelgv
        nelgt = abs(nelgs)
      endif
      call bcast(ldimr,ISIZE)
      call bcast(nelgs,ISIZE)
      call bcast(nelgv,ISIZE)
      call bcast(nelgt,ISIZE)
      ifre2 = .false.
      if (nelgs.lt.0) ifre2 = .true.

      call usrdat0

      if (nelgt.gt.350000 .and. .not.ifre2) 
     $   call exitti('Problem size requires .re2!$',1)

      if (ifre2) call read_re2_hdr(ifbswap, .true.) ! rank0 will open and read
      call chk_nel  ! make certain sufficient array sizes

      call mapelpr

      if (ifre2) then
        call read_re2_data(ifbswap, .true., .true., .true.)
      else
        maxrd = 32               ! max # procs to read at once
        mread = (np-1)/maxrd+1   ! mod param
        iread = 0                ! mod param
        x     = 0
        do i=0,np-1,maxrd
           call nekgsync()
           if (mod(nid,mread).eq.iread) then
              if (nid.ne.0) then
                open(UNIT=9,FILE=REAFLE,STATUS='OLD')
                call cscan(string,'MESH DATA',9)
                read(9,*) string
              endif 
              call rdmesh
              call rdcurve !  Curved side data
              call rdbdry  !  Boundary Conditions
              if (nid.ne.0) close(unit=9)
           endif
           iread = iread + 1
        enddo
      endif

C     Read Restart options / Initial Conditions / Drive Force
      CALL RDICDF
C     Read materials property data
      CALL RDMATP
C     Read history data
      CALL RDHIST
C     Read output specs
      CALL RDOUT
C     Read objects
      CALL RDOBJ

      call nekgsync()

C     End of input data, close read file.
      if(nid.eq.0) then
        close(unit=9)
        call echopar
        write(6,'(A,g13.5,A,/)')  ' done :: read .rea file ',
     $                             dnekclock()-etime0,' sec'
      endif

 99   call izero(boundaryID, size(boundaryID))
      call izero(boundaryIDt, size(boundaryIDt))

      ifld = 2 
      if(ifflow) ifld = 1
      do iel = 1,nelv
      do ifc = 1,2*ndim   
         boundaryID(ifc,iel) = bc(5,ifc,iel,ifld)
      enddo
      enddo

      ntmsh = 0
      do i=1,ldimt
         if(iftmsh(1+i)) ntmsh = ntmsh + 1 
      enddo

      if (ntmsh.gt.0) then
        do iel = 1,nelt
        do ifc = 1,2*ndim   
           boundaryIDt(ifc,iel) = bc(5,ifc,iel,2)
        enddo
        enddo
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine vrdsmsh
C
C=====================================================================
C     Verify that mesh and dssum are properly defined by performing
C        a direct stiffness operation on the X,Y and Z coordinates.
C     Note that periodic faces are not checked here.
C=====================================================================
C
      include 'SIZE'
      include 'TOTAL'
      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT),TB(LX1,LY1,LZ1,LELT)
     $           ,QMASK(LX1,LY1,LZ1,LELT),tmp(2)
      CHARACTER*3 CB

c      call  vrdsmshx  ! verify mesh topology

      if(nio.eq.0) write(*,*) 'verify mesh topology'

      IERR      = 0
      EPS       = 1.0e-04
      EPS       = 1.0e-03
      IFIELD    = 1
      if (nelgv.ne.nelgt .or. .not.ifflow) ifield = 2
      NXYZ1     = lx1*ly1*lz1
      NTOT      = lx1*ly1*lz1*NELT
      NFACES    = 2*ldim

      xmx = glmax(xm1,ntot)
      xmn = glmin(xm1,ntot)
      ymx = glmax(ym1,ntot)
      ymn = glmin(ym1,ntot)
      zmx = glmax(zm1,ntot)
      zmn = glmin(zm1,ntot)
      if (nio.eq.0) write(6,*) xmn,xmx,' Xrange'
      if (nio.eq.0) write(6,*) ymn,ymx,' Yrange'
      if (nio.eq.0) write(6,*) zmn,zmx,' Zrange'
c     return
C
C     First check - use 1/Multiplicity
C
      IF (IFHEAT) THEN 
         CALL COPY(TA,TMULT,NTOT)
      ELSE
         CALL COPY(TA,VMULT,NTOT)
      ENDIF
c
c     write(6,1) 
c    $(nid,'tab4',lglel(ie),(ta(k,1,1,ie),k=1,lx1*ly1),ie=1,nelt)
c   1 format(i3,a4,i3,16f5.2)
c
      CALL DSSUM(TA,lx1,ly1,lz1)
c
c     write(6,1) 
c    $(nid,'taaf',lglel(ie),(ta(k,1,1,ie),k=1,lx1*ly1),ie=1,nelt)
c
      CALL RONE (TB,NTOT)
      CALL SUB2 (TB,TA,NTOT)
      DO 1000 IE=1,NELT
      ieg=lglel(ie)
      DO 1000 IZ=1,lz1
      DO 1000 IY=1,ly1
      DO 1000 IX=1,lx1
         IF (ABS(TB(IX,IY,IZ,IE)).GT.EPS ) THEN
            WRITE(6,1005) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TA(IX,IY,IZ,IE),eps
c           WRITE(7,1005) IX,IY,IZ,IEG
c    $      ,XM1(IX,IY,IZ,IE),TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE)
c    $      ,QMASK(IX,IY,IZ,IE)
 1005       FORMAT(2X,'WARNING: DSSUM problem at:',/
     $            ,1X,'I,J,K,IE:',3I5,i12,/
     $            ,2X,'Near X =',3G16.8,', d:',2G16.8)
            IERR=4
         ENDIF
 1000 CONTINUE
C
C     Set up QMASK quickly to annihilate checks on periodic bc's
C
      CALL RONE(QMASK,NTOT)
      DO 100 IEL=1,NELT
      DO 100 IFACE=1,NFACES
         CB =CBC(IFACE,IEL,IFIELD)
         IF (CB.EQ.'P  '.or.cb.eq.'p  ') 
     $         CALL FACEV(QMASK,IEL,IFACE,0.0,lx1,ly1,lz1)
  100 CONTINUE
      CALL DSOP(QMASK,'MUL',lx1,ly1,lz1)

c      xxmin = glmin(xm1,ntot)
c      yymin = glmin(ym1,ntot)
c      zzmin = glmin(zm1,ntot)
c      xxmax = glmax(xm1,ntot)
c      yymax = glmax(ym1,ntot)
c      zzmax = glmax(zm1,ntot)
c      if (nio.eq.0) write(6,7) xxmin,yymin,zzmin,xxmax,yymax,zzmax
c    7 format('xyz minmx2:',6g13.5)



C
C     X-component
C
      CALL COPY(TA,XM1,NTOT)
      CALL COPY(TB,XM1,NTOT)
      CALL DSOP(TA,'MIN',lx1,ly1,lz1)
      CALL DSOP(TB,'MAX',lx1,ly1,lz1)
      CALL SUB2(TA,XM1,NTOT)
      CALL SUB2(TB,XM1,NTOT)
      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)
      DO 1100 IE=1,NELT
         XSCMAX = VLMAX(XM1(1,1,1,IE),NXYZ1)
         XSCMIN = VLMIN(XM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(XSCMAX-XSCMIN)
         SCAL2=ABS(XSCMAX)
         SCAL3=ABS(XSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         XSCALE = 1./SCAL1
         ieg=lglel(ie)
         DO 1100 IZ=1,lz1
         DO 1100 IY=1,ly1
         DO 1100 IX=1,lx1
         if (abs(ta(ix,iy,iz,ie)*xscale).gt.eps .or.
     $       abs(tb(ix,iy,iz,ie)*xscale).gt.eps ) then
            write(6,1105) ix,iy,iz,ieg
     $      ,xm1(ix,iy,iz,ie),ym1(ix,iy,iz,ie),zm1(ix,iy,iz,ie)
     $      ,tb(ix,iy,iz,ie),ta(ix,iy,iz,ie),XSCALE
 1105       format(1x,'WARNING1 Element mesh mismatch at:',/
     $            ,1x,'i,j,k,ie:',3i5,I12,/
     $            ,1X,'Near X =',3G16.8,', d:',3G16.8)
            ierr=1
         endif
 1100 CONTINUE
C
C     Y-component
C
      CALL COPY(TA,YM1,NTOT)
      CALL COPY(TB,YM1,NTOT)
      CALL DSOP(TA,'MIN',lx1,ly1,lz1)
      CALL DSOP(TB,'MAX',lx1,ly1,lz1)
      CALL SUB2(TA,YM1,NTOT)
      CALL SUB2(TB,YM1,NTOT)
      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)
      DO 1200 IE=1,NELT
         YSCMAX = VLMAX(YM1(1,1,1,IE),NXYZ1)
         YSCMIN = VLMIN(YM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(YSCMAX-YSCMIN)
         SCAL2=ABS(YSCMAX)
         SCAL3=ABS(YSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         YSCALE = 1./SCAL1
         ieg=lglel(ie)
         DO 1200 IZ=1,lz1
         DO 1200 IY=1,ly1
         DO 1200 IX=1,lx1
         IF (ABS(TA(IX,IY,IZ,IE)*YSCALE).GT.EPS .OR.
     $       ABS(TB(IX,IY,IZ,IE)*YSCALE).GT.EPS ) THEN
            WRITE(6,1205) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),yscale
 1205       FORMAT(1X,'WARNING2 Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',3I5,i12,/
     $            ,1X,'Near Y =',3G16.8,', d:',3G16.8)
            IERR=2
         ENDIF
 1200 CONTINUE
C
C     Z-component
C
      IF (IF3D) THEN
       CALL COPY(TA,ZM1,NTOT)
       CALL COPY(TB,ZM1,NTOT)
       CALL DSOP(TA,'MIN',lx1,ly1,lz1)
       CALL DSOP(TB,'MAX',lx1,ly1,lz1)
       CALL SUB2(TA,ZM1,NTOT)
       CALL SUB2(TB,ZM1,NTOT)
       CALL COL2(TA,QMASK,NTOT)
       CALL COL2(TB,QMASK,NTOT)
       DO 1300 IE=1,NELT
          ZSCMAX = VLMAX(ZM1(1,1,1,IE),NXYZ1)
          ZSCMIN = VLMIN(ZM1(1,1,1,IE),NXYZ1)
          SCAL1=ABS(ZSCMAX-ZSCMIN)
          SCAL2=ABS(ZSCMAX)
          SCAL3=ABS(ZSCMIN)
          SCAL1=MAX(SCAL1,SCAL2)
          SCAL1=MAX(SCAL1,SCAL3)
          ZSCALE = 1./SCAL1
          ieg=lglel(ie)
          DO 1300 IZ=1,lz1
          DO 1300 IY=1,ly1
          DO 1300 IX=1,lx1
          IF (ABS(TA(IX,IY,IZ,IE)*ZSCALE).GT.EPS .OR.
     $        ABS(TB(IX,IY,IZ,IE)*ZSCALE).GT.EPS ) THEN
           WRITE(6,1305) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),zscale
 1305       FORMAT(1X,'WARNING3 Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',3I5,i12,/
     $            ,1X,'Near Z =',3G16.8,', d:',3G16.8)
            IERR=3
          ENDIF
 1300  CONTINUE
      ENDIF

      ierr = iglsum(ierr,1)
      IF (IERR.gt.0) THEN
         if(nid.eq.0) WRITE(6,1400) 
 1400    FORMAT
     $   (' Mesh consistency check failed.  EXITING in VRDSMSH.')
            call exitt
      ENDIF
 
      tmp(1)=ierr
      CALL GOP(tmp,tmp(2),'M  ',1)
      IF (tmp(1).ge.4.0) THEN
         WRITE(6,1400) 
     $   (' Mesh consistency check failed.  EXITING in VRDSMSH.')
         call exitt
      ENDIF
 
      if(nio.eq.0) then
        write(6,*) 'done :: verify mesh topology'
        write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine vrdsmshx  ! verify mesh topology
C
C=====================================================================
C     Verify that mesh and dssum are properly defined by performing
C        a direct stiffness operation on the X,Y and Z coordinates.
C     Note that periodic faces are not checked here.
C=====================================================================
C
      include 'SIZE'
      include 'TOTAL'
      common /scrns/ tc(lx1,ly1,lz1,lelt),td(lx1,ly1,lz1,lelt)
     $             , ta(lx1,ly1,lz1,lelt),tb(lx1,ly1,lz1,lelt)
     $             , qmask(lx1,ly1,lz1,lelt)
      CHARACTER*3 CB
C
      IERR      = 0
      EPS       = 1.0e-04
      EPS       = 1.0e-03
      IFIELD    = 1
      IF (IFHEAT) IFIELD = 2
      NXYZ1     = lx1*ly1*lz1
      NTOT      = lx1*ly1*lz1*NELT
      NFACES    = 2*ldim

      xmx = glmax(xm1,ntot)
      xmn = glmin(xm1,ntot)
      ymx = glmax(ym1,ntot)
      ymn = glmin(ym1,ntot)
      zmx = glmax(zm1,ntot)
      zmn = glmin(zm1,ntot)
      if (nio.eq.0) write(6,*) xmn,xmx,' Xrange'
      if (nio.eq.0) write(6,*) ymn,ymx,' Yrange'
      if (nio.eq.0) write(6,*) zmn,zmx,' Zrange'
c     return
C
C     First check - use 1/Multiplicity
C
      IF (IFHEAT) THEN 
         CALL COPY(TA,TMULT,NTOT)
      ELSE
         CALL COPY(TA,VMULT,NTOT)
      ENDIF
c
c     write(6,1) 
c    $(nid,'tab4',lglel(ie),(ta(k,1,1,ie),k=1,lx1*ly1),ie=1,nelt)
c   1 format(i3,a4,i3,16f5.2)
c
      CALL DSSUM(TA,lx1,ly1,lz1)
c
c     write(6,1) 
c    $(nid,'taaf',lglel(ie),(ta(k,1,1,ie),k=1,lx1*ly1),ie=1,nelt)
c
      CALL RONE (TB,NTOT)
      CALL SUB2 (TB,TA,NTOT)
      DO 1000 IE=1,NELT
      ieg=lglel(ie)
      DO 1000 IZ=1,lz1
      DO 1000 IY=1,ly1
      DO 1000 IX=1,lx1
         IF (ABS(TB(IX,IY,IZ,IE)).GT.EPS ) THEN
            WRITE(6,1005) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TA(IX,IY,IZ,IE),eps
c           WRITE(7,1005) IX,IY,IZ,IEG
c    $      ,XM1(IX,IY,IZ,IE),TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE)
c    $      ,QMASK(IX,IY,IZ,IE)
 1005       FORMAT(2X,'WARNING: DSSUM problem at:',/
     $            ,1X,'I,J,K,IE:',3I5,i12,/
     $            ,2X,'Near X =',3G16.8,', d:',2G16.8)
            IERR=4
         ENDIF
 1000 CONTINUE
C
C     Set up QMASK quickly to annihilate checks on periodic bc's
C
      CALL RONE(QMASK,NTOT)
      DO 100 IEL=1,NELT
      DO 100 IFACE=1,NFACES
         CB =CBC(IFACE,IEL,IFIELD)
         IF (CB.EQ.'P  '.or.cb.eq.'p  ') 
     $         CALL FACEV(QMASK,IEL,IFACE,0.0,lx1,ly1,lz1)
  100 CONTINUE
      call dsop(QMASK,'MUL',lx1,ly1,lz1)

      xxmin = glmin(xm1,ntot)
      yymin = glmin(ym1,ntot)
      zzmin = glmin(zm1,ntot)
      xxmax = glmax(xm1,ntot)
      yymax = glmax(ym1,ntot)
      zzmax = glmax(zm1,ntot)
      if (nio.eq.0) write(6,7) xxmin,yymin,zzmin,xxmax,yymax,zzmax
    7 format('xyz minmx2:',6g13.5)


C
C     X-component
C
      call copy(ta,xm1,ntot)
      call copy(tb,xm1,ntot)
      call dsop(ta,'min',lx1,ly1,lz1)
      call dsop(tb,'max',lx1,ly1,lz1)

      call copy(tc,xm1,ntot)
      call copy(td,xm1,ntot)
      call dsop(tc,'min',lx1,ly1,lz1)
      call dsop(td,'max',lx1,ly1,lz1)

      xxmin = glmin(xm1,ntot)
      xxmax = glmax(xm1,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nio.eq.0) write(6,9) xxmin,yymin,zzmin,xxmax,yymax,zzmax
    9 format('xyz minmx3:',6g13.5)

      CALL SUB2(TA,XM1,NTOT)
      CALL SUB2(TB,XM1,NTOT)

      xxmin = glmin(qmask,ntot)
      xxmax = glmax(qmask,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nio.eq.0) write(6,19) xxmin,yymin,zzmin,xxmax,yymax,zzmax
   19 format('xyz minmx4:',6g13.5)

      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)

      xxmin = glmin(qmask,ntot)
      xxmax = glmax(qmask,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nio.eq.0) write(6,29) xxmin,yymin,zzmin,xxmax,yymax,zzmax
   29 format('xyz minmx5:',6g13.5)

      DO 1100 IE=1,NELT
         XSCMAX = VLMAX(XM1(1,1,1,IE),NXYZ1)
         XSCMIN = VLMIN(XM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(XSCMAX-XSCMIN)
         SCAL2=ABS(XSCMAX)
         SCAL3=ABS(XSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         XSCALE = 1./SCAL1
         ieg=lglel(ie)
         DO 1100 IZ=1,lz1
         DO 1100 IY=1,ly1
         DO 1100 IX=1,lx1
         if (abs(ta(ix,iy,iz,ie)*xscale).gt.eps .or.
     $       abs(tb(ix,iy,iz,ie)*xscale).gt.eps ) then
            write(6,1105) nid,ix,iy,iz,ie,ieg
     $      ,xm1(ix,iy,iz,ie),tc(ix,iy,iz,ie),td(ix,iy,iz,ie)
     $      ,ym1(ix,iy,iz,ie),zm1(ix,iy,iz,ie)
     $      ,ta(ix,iy,iz,ie),tb(ix,iy,iz,ie),xscale
     $      ,qmask(ix,iy,iz,ie)
 1105       format(i4.4,1x,'ie:',3i3,i10,i10,1p9e11.3)
c1105       format(i4.4,1x,'ie:',3i3,i6,1p9e11.3)
            ierr=1
            goto 1101
         endif
 1100 CONTINUE
 1101 CONTINUE

      xxmin = glmin(xm1,ntot)
      xxmax = glmax(xm1,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nio.eq.0) write(6,39) xxmin,yymin,zzmin,xxmax,yymax,zzmax
   39 format('xyz minmx5:',6g13.5)

c     ifvo = .true.
c     ifpo = .false.
c     ifto = .true.
c     call outpost(xm1,ta,tb,pr,qmask,'   ')
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine rotat2(xyz,angle,npts)
C
C     Rotate NPTS through ANGLE (in two directions IF3D).
C
      include 'SIZE'
      include 'INPUT'
      DIMENSION XYZ(3,1)
      COMMON /CTMP0/ RMTRX(3,3),RX(3,3),RZ(3,3),XYZN(3,10)
C
      SINA=SIN(ANGLE)
      COSA=COS(ANGLE)
      CALL RZERO(RX,9)
      CALL RZERO(RZ,9)
      RX(1,1)=COSA
      RX(2,2)=COSA
      RX(1,2)=SINA
      RX(2,1)=-SINA
      RX(3,3)=1.0
      IF (IF3D) THEN
         RZ(1,1)=COSA
         RZ(3,3)=COSA
         RZ(1,3)=SINA
         RZ(3,1)=-SINA
         RZ(2,2)=1.0
      ELSE
         RZ(1,1)=1.0
         RZ(2,2)=1.0
         RZ(3,3)=1.0
      ENDIF
      CALL MXM(RX,3,RZ,3,RMTRX,3)
C
C     Strip mine mxms in chunks of 10:
      DO 100 I=1,NPTS-10,10
         CALL MXM(RMTRX,3,XYZ(1,I),3,XYZN,10)
         CALL COPY(XYZ(1,I),XYZN,30)
  100 CONTINUE
      N10=MOD1(NPTS,10)
      I=NPTS-N10+1
      CALL RZERO(XYZN,30)
      IF (N10.GT.0) THEN
         CALL MXM(RMTRX,3,XYZ(1,I),3,XYZN,N10)
         CALL COPY(XYZ(1,I),XYZN,3*N10)
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine scale(xyzl,nl)
C
C     Rescale XYZL such that the mean value of IXX=IYY=IZZ for each element.
C
      include 'SIZE'
      include 'INPUT'
      DIMENSION XYZL(3,8,LELT)
      COMMON /CTMP0/ VO(LELT),XYZI(3,LELT),CG(3,LELT)
     $              ,TI(6),WORK(6)
C
C     Compute volumes -
C
      CALL VOLUME2(VO,XYZL,NL)
      VTOT=GLSUM (VO,NL)
C
C     Compute (weighted) average inertia for each element.
C
      NCRNR=2**ldim
      CALL RZERO(TI,6)
      DO 100 IL=1,NL
         VO0 = VO(IL)/VTOT
         CALL INRTIA(XYZI(1,IL),CG(1,IL),XYZL(1,1,IL),NCRNR,1)
         TI(1)=TI(1)+XYZI(1,IL)*VO0
         TI(2)=TI(2)+XYZI(2,IL)*VO0
         TI(3)=TI(3)+XYZI(3,IL)*VO0
         TI(4)=TI(4)+CG(1,IL)  *VO0
         TI(5)=TI(5)+CG(2,IL)  *VO0
         TI(6)=TI(6)+CG(3,IL)  *VO0
  100 CONTINUE
      CALL GOP(TI,WORK,'+  ',6)
      XI  =SQRT(TI(1))
      YI  =SQRT(TI(2))
      ZI  =1.0
      IF (IF3D) ZI=SQRT(TI(3))
C
C     Rescale ( & shift to a nearly mean zero )
C
      DO 200 IL=1,NL
      DO 200 IC=1,NCRNR
         XYZL(1,IC,IL)=(XYZL(1,IC,IL)-TI(4))/XI
         XYZL(2,IC,IL)=(XYZL(2,IC,IL)-TI(5))/YI
         XYZL(3,IC,IL)=(XYZL(3,IC,IL)-TI(6))/ZI
  200 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine inrtia(xyzi,cg,xyzl,n,itype)
C
C     Compute cg and inertia for a collection of unit point masses.
C     This is a global (multiprocessor) operation, only IF itype=2.
C
      DIMENSION XYZI(3),CG(3),XYZL(3,1)
      DIMENSION TI(4),WORK(4)
C
      TI(1)=0.0
      TI(2)=0.0
      TI(3)=0.0
      TI(4)=N
      DO 100 I=1,N
         TI(1)=TI(1)+XYZL(1,I)
         TI(2)=TI(2)+XYZL(2,I)
         TI(3)=TI(3)+XYZL(3,I)
  100 CONTINUE
      IF (ITYPE.EQ.2) CALL GOP(TI,WORK,'+  ',4)
      IF (TI(4).EQ.0.0) TI(4)=1.0
      CG(1)=TI(1)/TI(4)
      CG(2)=TI(2)/TI(4)
      CG(3)=TI(3)/TI(4)
C
      TI(1)=0.0
      TI(2)=0.0
      TI(3)=0.0
      DO 200 I=1,N
         TI(1)=TI(1)+( XYZL(1,I)-CG(1) )**2
         TI(2)=TI(2)+( XYZL(2,I)-CG(2) )**2
         TI(3)=TI(3)+( XYZL(3,I)-CG(3) )**2
  200 CONTINUE
      IF (ITYPE.EQ.2) CALL GOP(TI,WORK,'+  ',3)
      TI(1)=TI(1)/TI(4)
      TI(2)=TI(2)/TI(4)
      TI(3)=TI(3)/TI(4)
      IF (ITYPE.EQ.2) THEN
C        std. def'n of inertia.
         XYZI(1)=TI(2)+TI(3)
         XYZI(2)=TI(3)+TI(1)
         XYZI(3)=TI(1)+TI(2)
      ELSE
         XYZI(1)=TI(1)
         XYZI(2)=TI(2)
         XYZI(3)=TI(3)
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine volume2(vol,xyz,n)
      include 'SIZE'
      include 'INPUT'
      DIMENSION XYZ(3,2,2,2,1)
      DIMENSION VOL(1)
C
      DO 1000 IE=1,N
         VOL(IE)=0.0
         IF (IF3D) THEN        
           DO 20 K=1,2
           DO 20 J=1,2
           DO 20 I=1,2
              VOL1 = (XYZ(1,2,J,K,IE)-XYZ(1,1,J,K,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,I,J,2,IE)-XYZ(3,I,J,1,IE))
              VOL2 = (XYZ(1,2,J,K,IE)-XYZ(1,1,J,K,IE))
     $             * (XYZ(2,I,J,2,IE)-XYZ(2,I,J,1,IE))
     $             * (XYZ(3,I,2,K,IE)-XYZ(3,I,1,K,IE))
              VOL3 = (XYZ(1,I,2,K,IE)-XYZ(1,I,1,K,IE))
     $             * (XYZ(2,2,J,K,IE)-XYZ(2,1,J,K,IE))
     $             * (XYZ(3,I,J,2,IE)-XYZ(3,I,J,1,IE))
              VOL4 = (XYZ(1,I,J,2,IE)-XYZ(1,I,J,1,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,I,2,K,IE)-XYZ(3,I,1,K,IE))
              VOL5 = (XYZ(1,I,2,K,IE)-XYZ(1,I,1,K,IE))
     $             * (XYZ(2,I,J,2,IE)-XYZ(2,I,J,1,IE))
     $             * (XYZ(3,2,J,K,IE)-XYZ(3,1,J,K,IE))
              VOL6 = (XYZ(1,I,J,2,IE)-XYZ(1,I,J,1,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,2,J,K,IE)-XYZ(3,1,J,K,IE))
              VOL(IE) = VOL(IE)+VOL1+VOL2+VOL3+VOL4+VOL5+VOL6
   20      CONTINUE
           VOL(IE)=VOL(IE)/8.0
         ELSE
C     2-D:
            DO 40 J=1,2
            DO 40 I=1,2
              VOL1 = (XYZ(1,2,J,1,IE)-XYZ(1,1,J,1,IE))
     $             * (XYZ(2,I,2,1,IE)-XYZ(2,I,1,1,IE))
              VOL3 = (XYZ(1,I,2,1,IE)-XYZ(1,I,1,1,IE))
     $             * (XYZ(2,2,J,1,IE)-XYZ(2,1,J,1,IE))
              VOL(IE)=VOL(IE)+VOL1+VOL3
   40      CONTINUE
           VOL(IE)=VOL(IE)/4.0
         ENDIF
         VOL(IE)=ABS(VOL(IE))
 1000 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine findcg(cg,xyz,n)
C
C     Compute cg for N elements.
C
      include 'SIZE'
      DIMENSION CG(3,1),XYZ(3,8,1)
C
      NCRNR=2**ldim
      CALL RZERO(CG,3*N)
      DO 100 I =1,N
      DO 100 IC=1,NCRNR
         CG(1,I)=CG(1,I)+XYZ(1,IC,I)
         CG(2,I)=CG(2,I)+XYZ(2,IC,I)
         CG(3,I)=CG(3,I)+XYZ(3,IC,I)
  100 CONTINUE
      TMP=1.0/(NCRNR)
      CALL CMULT(CG,TMP,3*N)
      return
      end
c-----------------------------------------------------------------------
      subroutine divide(list1,list2,nl1,nl2,ifok,list,nl,xyzi,cg,WGT)
C
C     Divide the elements associated with this subdomain according to
C     the direction having the smallest moment of inertia (the "long"
C     direction).
C
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
C
      DIMENSION LIST(LELT),LIST1(LELT),LIST2(LELT)
      DIMENSION XYZI(3),CG(3,LELT),wgt(1)
      COMMON /CTMP0/ XCG(LELT),YCG(LELT),ZCG(LELT)
      REAL IXX,IYY,IZZ
      INTEGER WORK(2),WRK2(2)
      LOGICAL IFOK
C
C     Choose "long" direction:
C
      IXX=XYZI(1)
      IYY=XYZI(2)
      IZZ=XYZI(3)
      IF (IF3D) THEN
         IF (IXX.LE.IYY.AND.IXX.LE.IZZ) THEN
            DO 104 IE=1,NL
               XCG(IE)=CG(1,IE)
               YCG(IE)=CG(2,IE)
               ZCG(IE)=CG(3,IE)
  104       CONTINUE
         ELSEIF (IYY.LE.IXX.AND.IYY.LE.IZZ) THEN
            DO 106 IE=1,NL
               XCG(IE)=CG(2,IE)
               YCG(IE)=CG(3,IE)
               ZCG(IE)=CG(1,IE)
  106       CONTINUE
         ELSEIF (IZZ.LE.IXX.AND.IZZ.LE.IYY) THEN
            DO 108 IE=1,NL
               XCG(IE)=CG(3,IE)
               YCG(IE)=CG(1,IE)
               ZCG(IE)=CG(2,IE)
  108       CONTINUE
         ENDIF
      ELSE
C     2-D:
         IF (IXX.LE.IYY) THEN
            DO 114 IE=1,NL
               XCG(IE)=CG(1,IE)
               YCG(IE)=CG(2,IE)
  114       CONTINUE
         ELSE
            DO 116 IE=1,NL
               XCG(IE)=CG(2,IE)
               YCG(IE)=CG(1,IE)
  116       CONTINUE
         ENDIF
      ENDIF
      call col2(xcg,wgt,nl)
      call col2(ycg,wgt,nl)
      call col2(zcg,wgt,nl)
C
C     Find median value of CG to determine dividing point:
C
      XM=FMDIAN(XCG,NL,IFOK)
      YM=FMDIAN(YCG,NL,IFOK)
      ZM=0.0
      IF (IF3D) ZM=FMDIAN(ZCG,NL,IFOK)
C
C     Diagnostics
C
      IF (.NOT.IFOK) THEN
         WRITE(6,130) NID,NL,XM,YM,ZM
         DO 120 IL=1,NL
            WRITE(6,135) NID,IL,XCG(IL),YCG(IL),ZCG(IL)
  120    CONTINUE
  130    FORMAT(I10,'DIVIDE: NL,XM,YM,ZM',I3,3F12.5)
  135    FORMAT(I10,'DIVIDE: NID,IL,XC,YC,ZCG',I4,3F12.5)
      ENDIF
C
C=============================================================
C     Divide LIST into LIST1 (XCG < XM) and LIST2 (XCG>XM).
C=============================================================
C
      NL1=0
      NL2=0
      DO 200 IE=1,NL
         IF (XCG(IE).LT.XM) THEN
            NL1=NL1+1
            LIST1(NL1)=LIST(IE)
         ENDIF
         IF (XCG(IE).GT.XM) THEN
            NL2=NL2+1
            LIST2(NL2)=LIST(IE)
         ENDIF
         IF (XCG(IE).EQ.XM) THEN
C
C           We have to look at the other directions to arrive at
C           a unique subdivision algortithm.
C
C
C           More Diagnostics
C
            IF (.NOT.IFOK) WRITE(6,201) NID,IE,XCG(IE),XM
  201    FORMAT(I10,'DIVIDE: IE,XCG,XM:',I4,3F12.5)
C
            IF (YCG(IE).LT.YM) THEN
               NL1=NL1+1
               LIST1(NL1)=LIST(IE)
            ENDIF
            IF (YCG(IE).GT.YM) THEN
               NL2=NL2+1
               LIST2(NL2)=LIST(IE)
            ENDIF
            IF (YCG(IE).EQ.YM) THEN
C              look at 3rd direction.  
               IF (IF3D .AND. ZCG(IE).LT.ZM) THEN
                  NL1=NL1+1
                  LIST1(NL1)=LIST(IE)
               ELSE IF (IF3D .AND. ZCG(IE).GT.ZM) THEN
                  NL2=NL2+1
                  LIST2(NL2)=LIST(IE)
               ELSE 
C                 for 2- or 3-D intdeterminate case:
                  NL1=NL1+1
                  LIST1(NL1)=LIST(IE)
               ENDIF
            ENDIF
C
         ENDIF
  200 CONTINUE
C
C     Check for an even distribution (i.e. - not different by
C     more than 1):
C
      IFOK=.TRUE.  
      WORK(1)=NL1 
      WORK(2)=NL2
      CALL IGOP(WORK,WRK2,'+  ',2)
      IF (ABS(WORK(1)-WORK(2)).GT.1) IFOK=.FALSE.
C
      return
      end
c-----------------------------------------------------------------------
      subroutine bufchk(buf,n)
      integer n
      real buf(n)
      do i=1,n
         write(6,*) buf(i), ' whhhh'
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_xyz
      include 'SIZE'
      include 'TOTAL'
      integer e,f,eg

      do e=1,nelt
         eg = lglel(e)
         write(6,1) nid,eg,e,(cbc(f,e,1),f=1,6)
      enddo
    1 format(3i12,6(1x,a3),'  cbc')

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_nel
      include 'SIZE'
      include 'TOTAL'

      neltmx=np*lelt
      nelvmx=np*lelv

      neltmx=min(neltmx,lelg)
      nelvmx=min(nelvmx,lelg)

      nelgt = iglmax(nelgt,1)
      nelgv = iglmax(nelgv,1)

c     write(6,*) nid,' inside chk_nel',nelgt,neltmx,nelvmx

      if (nelgt.gt.neltmx.or.nelgv.gt.nelvmx) then
         if (nid.eq.0) then
          lelt_needed = nelgt/np
          if (mod(nelgt,np).ne.0) lelt_needed = lelt_needed + 1 
          write(6,12) lelt,lelg,lelt_needed,np,nelgt
 12         format(//,2X,'ABORT: Problem size too large!'
     $         ,/,2X
     $         ,/,2X,'This solver has been compiled for:'
     $         ,/,2X,'   number of elements/proc  (lelt):',i12
     $         ,/,2X,'   total number of elements (lelg):',i12
     $         ,/,2X
     $         ,/,2X,'Recompile with the following SIZE  parameters:'
     $         ,/,2X,'   lelt >= ',i12,'  for np = ',i12
     $         ,/,2X,'   lelg >= ',i12,/)
c           write(6,*)'help:',lp,np,nelvmx,nelgv,neltmx,nelgt
c           write(6,*)'help:',lelt,lelv,lelgv
         endif
         call exitt
      endif

      if(nelgt.gt.nelgt_max) then
        if(nid.eq.0) write(6,*)
     $               'ABORT: Total number of elements too large!',
     $               '       nel_max = ', nelgt_max 
        call exitt
      endif

      if (nelt.gt.lelt) then
        write(6,'(A,3I12)') 'ABORT: nelt>lelt!', nid, nelt, lelt
        call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cscan(sout,key,nk)

      character*132 sout,key
      character*132 string
      character*1  string1(132)
      equivalence (string1,string)
c
      do i=1,100000000
         call blank(string,132)
         read (nk,80,end=100,err=100) string
         call chcopy(sout,string,132)
c        write (6,*) string
         if (indx1(string,key,nk).ne.0) return
      enddo
  100 continue
c
   80 format(a132)
      return

      end
c-----------------------------------------------------------------------
c
c     NEW READER
c
      subroutine readat_big
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'CTIMER'

      logical ifbswap

      call read_re2_hdr(ifbswap, .true.)

      nelt = nelgt / np
      do i = 1, mod(nelgt, np)
        if (np-i.eq.nid) nelt = nelt + 1
      enddo
      if (nelt .gt. lelt) then
        call exitti('nelt > lelt!$',nelt)
      endif

      call chkParam

      call read_re2_data_big(ifbswap) ! Data mapped by mod(eg,np)

      call mapelpr_big       ! Test source code

      call nekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine read_re2_data_big(ifbswap) ! big .re2 reader
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'CTIMER'

      logical ifbswap

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
 

      etime0 = dnekclock_sync()

                  ibc = 2
      if (ifflow) ibc = 1

                  nfldt = 1
      if (ifheat) nfldt = 2+npscal
      if (ifmhd ) nfldt = 2+npscal+1

      ! first field to read
      if (param(33).gt.0) ibc = int(param(33))

      ! number of fields to read
      if (param(32).gt.0) then
        nfldt = ibc + int(param(32)) - 1
        nfldt = max(nfldt,1) 
        if (nelgt.gt.nelgv) nfldt = max(nfldt,2) 
      endif

      call blank(cbc,3*size(cbc))
      call rzero(bc ,size(bc))

#ifndef NOMPIIO
      call fgslib_crystal_setup(cr_re2,nekcomm,np)

      call byte_open_mpi(re2fle,fh_re2,.true.,ierr)
      call err_chk(ierr,' Cannot open .re2 file!$')

      call reado_re2_mesh(ifbswap)
      call reado_re2_curve(ifbswap)
      do ifield = ibc,nfldt
        call reado_re2_bc(cbc(1,1,ifield),bc(1,1,1,ifield),ifbswap)
      enddo

      call fgslib_crystal_free(cr_re2)
      call byte_close_mpi(fh_re2,ierr)
#else

      call exitti('No serial support for big mesh read! P=$',np)

#endif

      etime_t = dnekclock_sync() - etime0
      if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: read .re2 file   ',
     &                   etime_t, ' sec'

      return
      end
c-----------------------------------------------------------------------
      subroutine reado_re2_mesh(ifbswap) ! version 3 of .re2 reader

c     "read only" --- redistribute via q = mod(eg,np)

      include 'SIZE'
      include 'TOTAL'

      logical ifbswap

      parameter(nrmax = lelt)             ! maximum number of records
      parameter(lrs   = 1+ldim*(2**ldim)) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer e,eg,ind(nrmax)

      integer         bufr(li-2,nrmax)
      common /scrns/  bufr

      integer         vi  (li  ,nrmax)
      common /ctmp1/  vi

      integer*8       lre2off_b,dtmp8
      integer*8       nrg      

      nrg       = nelgt
      nr        = nelt
      irankoff  = igl_running_sum(nr) - nr
      dtmp8     = irankoff
      re2off_b  = 84 ! set initial offset (hdr + endian)
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      ! read coordinates from file
      nwds4r = nr*lrs4
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(bufr,nwds4r,-1,fh_re2,ierr)
      re2off_b = re2off_b + nrg*4*lrs4
      if (ierr.gt.0) goto 100

      if (nio.eq.0) write(6,*) 'reading mesh '

      ! pack buffer
      do i = 1,nr
         jj      = (i-1)*lrs4 + 1
         eg      = irankoff + i ! elements are stored in global order
         vi(1,i) = mod(eg,np)
         vi(2,i) = eg
         call icopy(vi(3,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      n   = nr
      key = 1 
      call fgslib_crystal_tuple_transfer(cr_re2,n,nrmax,vi,li,
     &   vl,0,vr,0,key)

      ! unpack buffer
      ierr = 0
      if (n.gt.nrmax) then
         ierr = 1
         goto 100
      endif

      do i = 1,n
         lglel(i)=vi(2,i)  ! List of global element numbers in v(2,:)
      enddo
      call isort(lglel,ind,n)

      do e = 1,n
         i = ind(e)
         call icopy     (bufr,vi(3,i),lrs4)
         call buf_to_xyz(bufr,e,ifbswap,ierr)
      enddo
      nelt = n

 100  call err_chk(ierr,'Error reading .re2 mesh$')

      return
      end
c-----------------------------------------------------------------------
      subroutine reado_re2_curve(ifbswap)

      include 'SIZE'
      include 'TOTAL'

      logical ifbswap

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      parameter(nrmax = 12*lelt) ! maximum number of records
      parameter(lrs   = 2+1+5)   ! record size: eg iside curve(5) ccurve
      parameter(li    = 2*lrs+1)

      integer         bufr(li-1,nrmax)
      common /scrns/  bufr

      integer         vi  (li  ,nrmax)
      common /ctmp1/  vi

      integer*8       lre2off_b,dtmp8
      integer*8       nrg,nr
      integer*4       nrg4(2)
     
      integer*8       i8gl_running_sum

      integer e,eg,esort(nrmax),ind(nrmax)

      ! read total number of records
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(nrg4,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      if(wdsizi.eq.8) then
         if(ifbswap) call byte_reverse8(nrg4,nwds4r,ierr)
         call copy(dnrg,nrg4,1)
         nrg = dnrg
      else
         if(ifbswap) call byte_reverse (nrg4,nwds4r,ierr)
         nrg = nrg4(1)
      endif
      re2off_b = re2off_b + 4*nwds4r

      if(nrg.eq.0) return

      ! read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      do i = 0,mod(nrg,dtmp8)-1
         if(i.eq.nid) nr = nr + 1
      enddo
      dtmp8     = i8gl_running_sum(int(nr,8)) - nr
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      re2off_b = re2off_b + nrg*4*lrs4

      if(nio.eq.0) write(6,'(A,I20)') ' reading curved sides   ', nrg

      nwds4r = nr*lrs4
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(bufr,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      ! pack buffer
      do i = 1,nr
         jj = (i-1)*lrs4 + 1

         if(ifbswap) then 
           lrs4s = lrs4 - wdsizi/4 ! words to swap (last is char)
           if(wdsizi.eq.8) call byte_reverse8(bufr(jj,1),lrs4s,ierr)
           if(wdsizi.eq.4) call byte_reverse (bufr(jj,1),lrs4s,ierr)
         endif

         eg = bufr(jj,1)
         if (wdsizi.eq.8) call copyi4(eg,bufr(jj,1),1)

         if(eg.le.0 .or. eg.gt.nelgt) goto 100
         vi(1,i) = mod(eg,np)

         call icopy (vi(2,i),bufr(jj,1),lrs4)  ! First entry is "eg"
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      n    = nr
      key  = 1
      call fgslib_crystal_tuple_transfer(cr_re2,n,nrmax,vi,li,vl,0,vr,0,
     &                                   key)

      if (wdsizi.eq.8) then
        do i = 1,n
          call copyi4(esort(i),vi(2,i),1)
        enddo
      else if (wdsizi.eq.4) then
        do i = 1,n
          esort(i) = vi(2,i)
        enddo
      endif
      call isort(esort,ind,n)

      if(n.gt.nrmax) goto 100

      ! unpack buffer
      do k = 1,n
         i  = ind(k)
         call find_lglel_ind(e, esort(k), nelt)
         ! Face f is dereferenced from vi(2,:)
         call buf_to_curve_loc(e,vi(2,i))
      enddo

      return

 100  ierr = 1
      call err_chk(ierr,'Error reading .re2 curved data$')

      end
c-----------------------------------------------------------------------
      subroutine read_re2_curve_v2(nvi, vi, ifbswap)
      include 'SIZE'
      include 'TOTAL'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      parameter(nrmax = 12*lelt) ! maximum number of records
      parameter(lrs   = 8)   ! record size: eg iside curve(5) ccurve
      parameter(li    = 2*lrs+1)

      integer nvi, vi(li, nrmax)
      logical ifbswap

      integer         bufr(li-1,nrmax)
      common /scrns/  bufr

      integer*8       lre2off_b,dtmp8
      integer*8       nrg,nr
      integer*4       nrg4(2)
     
      integer*8       i8gl_running_sum

      integer e, eg, ierr

      ! read total number of records
      nvi = 0
      ierr = 0
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(nrg4,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      if(wdsizi.eq.8) then
         if(ifbswap) call byte_reverse8(nrg4,nwds4r,ierr)
         call copy(dnrg,nrg4,1)
         nrg = dnrg
      else
         if(ifbswap) call byte_reverse (nrg4,nwds4r,ierr)
         nrg = nrg4(1)
      endif
      re2off_b = re2off_b + 4*nwds4r

      if(nrg.eq.0) return

      ! read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      do i = 0,mod(nrg,dtmp8)-1
         if(i.eq.nid) nr = nr + 1
      enddo
      dtmp8     = i8gl_running_sum(int(nr,8)) - nr
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      re2off_b = re2off_b + nrg*4*lrs4

      if(nid.eq.0) write(6,'(A,I20)') ' reading curved sides   ', nrg

      nwds4r = nr*lrs4
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(bufr,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      ! pack buffer
      do i = 1, nr
         jj = (i-1)*lrs4 + 1

         if(ifbswap) then
           lrs4s = lrs4 - wdsizi/4 ! words to swap (last is char)
           if(wdsizi.eq.8) call byte_reverse8(bufr(jj,1),lrs4s,ierr)
           if(wdsizi.eq.4) call byte_reverse (bufr(jj,1),lrs4s,ierr)
         endif

         eg = bufr(jj,1)
         if (wdsizi.eq.8) call copyi4(eg,bufr(jj,1),1)
         if(eg.le.0 .or. eg.gt.nelgt) goto 100

         call get_dest_proc(vi(1,i), eg, np, nelgt)
         ! First entry is "eg"
         call icopy (vi(2,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nvi = nr
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nvi, nrmax, vi, li,
     $  vl, 0, vr, 0, key)

      if(nvi.gt.nrmax) then
        ierr = 1
      endif

 100  call err_chk(ierr,'Error reading .re2 curved data$')
      return
      end
c-----------------------------------------------------------------------
      subroutine reado_re2_bc(cbl,bl,ifbswap)

      include 'SIZE'
      include 'TOTAL'

      character*3  cbl(  6,lelt)
      real         bl (5,6,lelt)
      logical      ifbswap

      parameter(nrmax = 6*lelt) ! maximum number of records
      parameter(lrs   = 2+1+5)  ! record size: eg iside bl(5) cbl
      parameter(li    = 2*lrs+1)

      integer         bufr(li-1,nrmax)
      common /scrns/  bufr

      integer         vi  (li  ,nrmax)
      common /ctmp1/  vi

      integer*8       lre2off_b,dtmp8
      integer*8       nrg,nr
      integer*4       nrg4(2)
      integer*8       i8gl_running_sum 

      integer e,eg,esort(nrmax),ind(nrmax)

      ! read total number of records
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(nrg4,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      if(wdsizi.eq.8) then
         if(ifbswap) call byte_reverse8(nrg4,nwds4r,ierr)
         call copy(dnrg,nrg4,1)
         nrg = dnrg
      else
         if(ifbswap) call byte_reverse (nrg4,nwds4r,ierr)
         nrg = nrg4(1)
      endif
      re2off_b = re2off_b + 4*nwds4r

      if(nrg.eq.0) return

      ! read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      do i = 0,mod(nrg,dtmp8)-1
         if(i.eq.nid) nr = nr + 1
      enddo
      dtmp8     = i8gl_running_sum(int(nr,8)) - nr
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      re2off_b = re2off_b + nrg*4*lrs4


      if(nio.eq.0) write(6,'(A,I20,A,I3)') 
     $             ' reading boundary faces ', nrg, 
     $             ' for ifield ', ifield

      nwds4r = nr*lrs4
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(bufr,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      ! pack buffer
      do i = 1,nr
         jj = (i-1)*lrs4 + 1

         if(ifbswap) then 
           lrs4s = lrs4 - wdsizi/4 ! words to swap (last is char)
           if(wdsizi.eq.8) call byte_reverse8(bufr(jj,1),lrs4s,ierr)
           if(wdsizi.eq.4) call byte_reverse (bufr(jj,1),lrs4s,ierr)
         endif

         eg = bufr(jj,1)
         if(wdsizi.eq.8) call copyi4(eg,bufr(jj,1),1)

         if(eg.le.0 .or. eg.gt.nelgt) goto 100
         vi(1,i) = mod(eg,np)

         call icopy (vi(2,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      n    = nr
      key  = 1

      call fgslib_crystal_tuple_transfer(cr_re2,n,nrmax,vi,li,vl,0,vr,0,
     &                                   key)

      if (wdsizi.eq.8) then
        do i = 1,n
          call copyi4(esort(i),vi(2,i),1)
        enddo
      else if (wdsizi.eq.4) then
        do i = 1,n
          esort(i) = vi(2,i)
        enddo
      endif
      call isort(esort,ind,n)

      ! fill up with default
      do e=1,nelt
      do k=1,6
         cbl(k,e) = 'E  '
      enddo
      enddo

      ! unpack buffer
      if (n.gt.nrmax) goto 100

      do k = 1,n
         i  = ind(k)
         call find_lglel_ind(e, esort(k), nelt)
         call buf_to_bc_loc(e,cbl,bl,vi(2,i))
      enddo

      return

 100  ierr = 1
      call err_chk(ierr,'Error reading .re2 boundary data$')

      end
c-----------------------------------------------------------------------
      subroutine read_re2_bc_v2(nvi, vi, cbl, bl, ifbswap)
      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = 6*lelt)
      parameter(lrs = 8)
      parameter(li = 2*lrs + 1)

      integer nvi, vi(li, nrmax)
      character*3 cbl(6, lelt)
      real bl(5, 6, lelt)
      logical ifbswap

      integer         bufr(li-1,nrmax)
      common /scrns/  bufr

      integer*8       lre2off_b,dtmp8
      integer*8       nrg,nr
      integer*4       nrg4(2)
      integer*8       i8gl_running_sum

      integer i, e, k, eg, ierr, key

      ! read total number of records
      nvi       = 0
      ierr      = 0
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(nrg4,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      if(wdsizi.eq.8) then
         if(ifbswap) call byte_reverse8(nrg4,nwds4r,ierr)
         call copy(dnrg,nrg4,1)
         nrg = dnrg
      else
         if(ifbswap) call byte_reverse (nrg4,nwds4r,ierr)
         nrg = nrg4(1)
      endif
      re2off_b = re2off_b + 4*nwds4r

      if(nrg.eq.0) return

      ! read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      do i = 0,mod(nrg,dtmp8)-1
         if(i.eq.nid) nr = nr + 1
      enddo
      dtmp8     = i8gl_running_sum(int(nr,8)) - nr
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      re2off_b = re2off_b + nrg*4*lrs4


      if(nio.eq.0) write(6,'(A,I20,A,I3)')
     $             ' reading boundary faces ', nrg,
     $             ' for ifield ', ifield

      nwds4r = nr*lrs4
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(bufr,nwds4r,-1,fh_re2,ierr)
      if(ierr.gt.0) goto 100

      ! pack buffer
      do i = 1, nr
         jj = (i-1)*lrs4 + 1

         if(ifbswap) then
           lrs4s = lrs4 - wdsizi/4 ! words to swap (last is char)
           if(wdsizi.eq.8) call byte_reverse8(bufr(jj,1),lrs4s,ierr)
           if(wdsizi.eq.4) call byte_reverse (bufr(jj,1),lrs4s,ierr)
         endif

         eg = bufr(jj,1)
         if(wdsizi.eq.8) call copyi4(eg,bufr(jj,1),1)
         if(eg.le.0 .or. eg.gt.nelgt) then
           ierr = 1
           goto 100
         endif

         call get_dest_proc(vi(1,i), eg, np, nelgt)
         call icopy (vi(2,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nvi = nr
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nvi, nrmax, vi, li,
     $  vl, 0, vr, 0, key)

      if (nvi.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      ! fill up with default
      do e=1,nelt
      do k=1,6
         cbl(k,e) = 'E  '
      enddo
      enddo

 100  call err_chk(ierr, 'Error reading .re2 boundary data$')
      return
      end
c-----------------------------------------------------------------------
      subroutine find_lglel_ind(e, eg, nl)
        include 'SIZE'
        include 'TOTAL'

        integer e, eg, nl
        call find_lglel_ind_binary(e, eg, nl)
        return
      end
c-----------------------------------------------------------------------
      subroutine find_lglel_ind_linear(e, eg, nl)
        include 'SIZE'
        include 'TOTAL'

        integer e, eg, nl, i, ierr

        ierr = 1
        do i = 1, nl
          if (lglel(i).eq.eg) then
            ierr = 0
            e = i
            return
          endif
        enddo

        call err_chk(ierr,'Error finding index in lglel$')
      end
c-----------------------------------------------------------------------
      subroutine find_lglel_ind_binary(e, eg, nl)
        include 'SIZE'
        include 'TOTAL'

        integer e, eg, nl
        integer i, j, mid, ierr

        ierr = 1
        if (nl.eq.0) then
          goto 100
        endif
        if ((eg.lt.lglel(1)).or.(eg.gt.lglel(nl))) then
          goto 100
        endif

        i = 1
        j = nl
        do while (i.lt.j)
          mid = (i + j)/2
          if (eg.gt.lglel(mid)) then
            i = mid + 1
          else
            j = mid
          endif
          if (eg.eq.lglel(mid)) then
            e = mid
            return
          endif
        enddo

        if (eg.eq.lglel(i)) then
          e = i
          return
        endif
        if (eg.eq.lglel(j)) then
          e = j
          return
        endif

 100    call err_chk(ierr, 'Error finding index in lglel binary$')
        return
      end
c-----------------------------------------------------------------------
      subroutine find_sorted_ind(e, eg, nl, sorted)
        include 'SIZE'
        include 'TOTAL'

        integer e, eg, nl, sorted(1)
        call find_sorted_ind_binary(e, eg, nl, sorted)
        return
      end
c-----------------------------------------------------------------------
c  Double check this
      subroutine find_sorted_ind_binary(e, eg, nl, sorted)
        include 'SIZE'
        include 'TOTAL'

        integer e, eg, nl, sorted(lelt)
        integer i, j, mid, ierr

        ierr = 1
        if (nl.eq.0) then
          goto 100
        endif
        if ((eg.lt.sorted(1)).or.(eg.gt.sorted(nl))) then
          goto 100
        endif

        i = 1
        j = nl
        do while (i.lt.j)
          mid = (i + j)/2
          if (eg.gt.sorted(mid)) then
            i = mid + 1
          else
            j = mid
          endif
          if (eg.eq.sorted(mid)) then
            e = mid
            return
          endif
        enddo

        if (eg.eq.sorted(i)) then
          e = i
          return
        endif
        if (eg.eq.sorted(j)) then
          e = j
          return
        endif

 100    call err_chk(ierr, 'Error finding index in sorted$')
        return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_curve_loc(e,buf)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'

      integer e,eg,f,buf(30)

      if(wdsizi.eq.8) then
        call copyi4(eg,buf(1),1) !1-2

        call copyi4(f,buf(3),1) !3-4

        call copy  ( curve(1,f,e),buf(5) ,5) !5--14
        call chcopy(ccurve(  f,e),buf(15),1)!15
      else
        eg = buf(1)
        f  = buf(2)

        call copy4r( curve(1,f,e),buf(3),5)
        call chcopy(ccurve(f,e)  ,buf(8),1)
      endif

c     write(6,1) eg,e,f,(curve(k,f,e),k=1,5),ccurve(f,e)
c   1 format(2i7,i3,5f10.3,1x,a1,'ccurve')

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_bc_loc(e,cbl,bl,buf)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'

      character*3 cbl(6,lelt)
      real         bl(5,6,lelt)

      integer e,eg,f,buf(30)

      if(wdsizi.eq.8) then
        call copyi4(eg,buf(1),1) !1-2

        call copyi4(f,buf(3),1) !3-4

        call copy  (bl(1,f,e),buf(5),5) !5--14
        call chcopy(cbl( f,e),buf(15),3)!15-16

        if(nelgt.ge.1000000.and.cbl(f,e).eq.'P  ')
     $   call copyi4(bl(1,f,e),buf(5),1) !Integer assign connecting P element

      else
        eg = buf(1)
        f  = buf(2)

        call copy4r ( bl(1,f,e),buf(3),5)
        call chcopy (cbl(  f,e),buf(8),3)

        if (nelgt.ge.1000000.and.cbl(f,e).eq.'P  ')
     $     bl(1,f,e) = buf(3) ! Integer assign of connecting periodic element
      endif

c      write(6,1) eg,e,f,cbl(f,e),' CBC',nid
c  1   format(2i8,i4,2x,a3,a4,i8)

      return
      end
c-----------------------------------------------------------------------
      subroutine transfer_vertices(vertex, loc_to_glob_nid)

      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = lelt)    ! maximum number of records
      parameter(lrs   = 2**ldim) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer*8 vertex(2**ldim, lelt)
      integer loc_to_glob_nid(lelt)

      integer         bufr(li - 2, nrmax)
      common /scrns/  bufr

      integer         vi  (li    , nrmax)
      common /ctmp1/  vi

      integer e, eg, ind(nrmax), nr, key

      lrs4      = lrs*wdsizi/4

      do e = 1, nelt
        vi(1, e) = loc_to_glob_nid(e)
        vi(2, e) = lglel(e)
        call vtx_to_buf(vi(3, e), vertex(1, e))
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nr = nelt
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     &  vl, 0, vr, 0, key)

      ! unpack buffer
      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      ! List of global element numbers in v(2,:)
      nelt = 0
      nelv = 0
      do i = 1, nr
        lglel(i) = vi(2, i)
        if (lglel(i).le.nelgv) then
          nelv = nelv + 1
        else
          nelt = nelt + 1
        endif
      enddo
      nelt = nelv + nelt
      call isort(lglel, ind, nr)

      do e = 1, nr
         i = ind(e)
         call buf_to_vtx(vertex(1, e), vi(3, i))
      enddo

 100  call err_chk(ierr, 'Error transferring vertices$')
      return
      end
c-----------------------------------------------------------------------
      subroutine transfer_vertices_v2(vertex, loc_to_glo_nid)
      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = lelt)    ! maximum number of records
      parameter(lrs   = 2**ldim) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer*8 vertex(2**ldim, lelt)
      integer loc_to_glo_nid(lelt)

      integer         bufr(li - 2, nrmax)
      common /scrns/  bufr

      integer         vi  (li    , nrmax)
      common /ctmp1/  vi

      integer e, eg, ind(nrmax), nr, key

      lrs4 = lrs*wdsizi/4

      do e = 1, nelt
        vi(1, e) = loc_to_glo_nid(e)
        vi(2, e) = lglel(e)
        call vtx_to_buf(vi(3, e), vertex(1, e))
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nr = nelt
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     &  vl, 0, vr, 0, key)

      ! unpack buffer
      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      ! List of global element numbers in v(2,:)
      nelt = 0
      nelv = 0
      do i = 1, nr
        lglel(i) = vi(2, i)
        if (lglel(i).le.nelgv) then
          nelv = nelv + 1
        else
          nelt = nelt + 1
        endif
      enddo
      nelt = nelv + nelt
      call isort(lglel, ind, nelt)

      do e = 1, nelt
         i = ind(e)
         call buf_to_vtx(vertex(1, e), vi(3, i))
      enddo

 100  call err_chk(ierr, 'Error transferring vertices$')
      return
      end
c-----------------------------------------------------------------------
      subroutine vtx_to_buf(buf, vtx)
        include 'SIZE'
        include 'TOTAL'

        integer*8 buf(2**ldim), vtx(2**ldim)
        integer e

        do e = 1, 2**ndim
          buf(e) = vtx(e)
        enddo

        return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_vtx(vtx, buf)
        include 'SIZE'
        include 'TOTAL'

        integer*8 vtx(2**ldim), buf(2**ldim)
        integer i

        do i = 1, 2**ldim
          vtx(i) = buf(i)
        enddo

        return
      end
c-----------------------------------------------------------------------
      subroutine transfer_con(wk, nwk, nelti, nlvi, ierr)

      include 'SIZE'
      include 'TOTAL'

      integer nwk, nelti, nlvi, ierr, wk(nwk)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      parameter(nrmax = lelt)    ! maximum number of records
      parameter(lrs   = 2**ldim) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer         vi(li, nrmax)
      common /ctmp1/  vi

      integer e, i, cnt, ind(nrmax), eg(nrmax), nr, key

      i = 1
      do e = 1, nelti
        vi(1, e) = mod(wk(i), np)
        vi(2, e) = wk(i)
        call vtx_to_buf_4(vi(3, e), wk(i + 1), nlvi)
        i = i + nlvi + 1
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nr = nelti
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     &  vl, 0, vr, 0, key)

      ! unpack buffer
      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      ! List of global element numbers in v(2,:)
      do i = 1, nr
        eg(i) = vi(2, i)
      enddo
      call isort(eg, ind, nr)

      cnt = 1
      do e = 1, nr
        wk(cnt) = eg(e)
        i = ind(e)
        call buf_to_vtx_4(wk(cnt + 1), vi(3, i), nlvi)
        cnt = cnt + nlvi + 1
      enddo

 100  call err_chk(ierr, 'Error transferring vertices$')
      return
      end
c-----------------------------------------------------------------------
      subroutine vtx_to_buf_4(buf, vtx, nlvi)
        include 'SIZE'
        include 'TOTAL'

        integer buf(2**ldim), vtx(2**ldim), nlvi
        integer e

        do e = 1, nlvi
          buf(e) = vtx(e)
        enddo

        return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_vtx_4(vtx, buf, nlvi)
        include 'SIZE'
        include 'TOTAL'

        integer vtx(2**ldim), buf(2**ldim), nlvi
        integer i

        do i = 1, nlvi
          vtx(i) = buf(i)
        enddo

        return
      end
c-----------------------------------------------------------------------
      subroutine transfer_re2_mesh(loc_to_glob_nid, lglelo, nelto)

      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = lelt)             ! maximum number of records
      parameter(lrs   = 1+ldim*(2**ldim)) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer loc_to_glob_nid(lelt), lglelo(lelt), nelto
      integer i, e, sorted(nrmax), ind(nrmax), nr, key

      integer         bufr(li - 2, nrmax)
      common /scrns/  bufr

      integer         vi  (li    , nrmax)
      common /ctmp1/  vi

      lrs4      = lrs*wdsizi/4

      do e = 1, nelto
        vi(1, e) = loc_to_glob_nid(e)
        vi(2, e) = lglelo(e)
        call xyz_to_buf(vi(3, e), e)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nr = nelto
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     &  vl, 0, vr, 0, key)

      ! unpack buffer
      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      ! List of global element numbers in v(2,:)
      do i = 1, nr
         sorted(i) = vi(2, i)
      enddo
      call isort(sorted, ind, nr)

      do e = 1, nr
         i = ind(e)
         call icopy(bufr, vi(3, i), lrs4)
         call buf_to_xyz(bufr, e, .false., ierr)
      enddo

 100  call err_chk(ierr, 'Error transferring .re2 mesh$')
      return
      end
c-----------------------------------------------------------------------
      subroutine transfer_re2_mesh_v2(loc_to_glo_nid, lglelo, nelto)
      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = lelt)             ! maximum number of records
      parameter(lrs   = 1+ldim*(2**ldim)) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer loc_to_glo_nid(lelt), lglelo(lelt), nelto
      integer e, sorted(nrmax), ind(nrmax), nr, key

      integer         bufr(li - 2, nrmax)
      common /scrns/  bufr

      integer         vi  (li    , nrmax)
      common /ctmp1/  vi

      lrs4 = lrs*wdsizi/4

      do e = 1, nelto
        vi(1, e) = loc_to_glo_nid(e)
        vi(2, e) = lglelo(e)
        call xyz_to_buf(vi(3, e), e)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      nr = nelto
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     &  vl, 0, vr, 0, key)

      ! unpack buffer
      ierr = 0
      if (nr.gt.nrmax .or. nr.ne.nelt) then
        ierr = 1
        goto 100
      endif

      ! List of global element numbers in v(2,:)
      do e = 1, nr
         sorted(e) = vi(2, e)
      enddo
      call isort(sorted, ind, nr)

      do e = 1, nr
         call icopy(bufr, vi(3, ind(e)), lrs4)
         call buf_to_xyz(bufr, e, .false., ierr)
      enddo

 100  call err_chk(ierr, 'Error transferring .re2 mesh$')
      return
      end
c-----------------------------------------------------------------------
      subroutine xyz_to_buf(buf, e)
        include 'SIZE'
        include 'TOTAL'

        integer buf(0:49), e

        if (wdsizi.eq.8) then
          call icopy48(buf(0), igroup(e), 1)
          if (ndim.eq.3) then
            call copy(buf( 2), xc(1, e), 8)
            call copy(buf(18), yc(1, e), 8)
            call copy(buf(34), zc(1, e), 8)
          else
            call copy(buf( 2), xc(1, e), 4)
            call copy(buf(10), yc(1, e), 4)
          endif
        else
          buf(0) = igroup(e)
          if (ndim.eq.3) then
            call copyX4(buf( 1), xc(1, e), 8)
            call copyX4(buf( 9), yc(1, e), 8)
            call copyX4(buf(17), zc(1, e), 8)
          else
            call copyX4(buf( 1), xc(1, e), 4)
            call copyX4(buf( 5), yc(1, e), 4)
          endif
        endif

        return
      end
c-----------------------------------------------------------------------
      subroutine transfer_re2_bc(cbl, bl, loc_to_glob_nid, lglelo,
     $  nelto)

      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = 6*lelt)
      parameter(lrs = 8)
      parameter(li = 2*lrs + 1)

      character*3 cbl(2*ldim, lelt)
      real bl(5, 2*ldim, lelt)
      integer loc_to_glob_nid(lelt), lglelo(lelt), nelto

      integer e, eg, f, cnt, nr, key, ierr, nf

      integer         vi(li, nrmax)
      common /ctmp1/  vi

      cnt = 0
      nf = 2*ndim
      do e = 1, nelto
        eg = lglelo(e)
        do f = 1, nf
          cnt = cnt + 1
          vi(1, cnt) = loc_to_glob_nid(e)
          call bc_to_buf(vi(2, cnt), cbl(f, e), bl(1, f, e), eg, f)
        enddo
      enddo

      nr = cnt
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     $  vl, 0, vr, 0, key)

      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      do i = 1,nr
         call buf_to_bc_big(cbl, bl, vi(2, i))
      enddo

 100  call err_chk(ierr, 'Error transferring .re2 boundary data$')
      return
      end
c-----------------------------------------------------------------------
      subroutine transfer_re2_bc_v2(nvi, vi, cbl, bl, loc_to_glo_nid,
     $  lglelo, nelto)
      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = 6*lelt)
      parameter(lrs = 8)
      parameter(li = 2*lrs + 1)

      integer nvi, vi(li, nrmax)
      character*3 cbl(6, lelt)
      real bl(5, 6, lelt)
      integer loc_to_glo_nid(lelt), lglelo(lelt), nelto

      integer i, e, nr, key, ierr
      integer eg(nrmax)

      ierr = 0
      if (wdsizi.eq.8) then
        do i = 1, nvi
          call copyi4(eg(i), vi(2, i), 1)
        enddo
      else if (wdsizi.eq.4) then
        do i = 1, nvi
          eg(i) = vi(2, i)
        enddo
      endif

      do i = 1, nvi
        call find_sorted_ind(e, eg(i), nelto, lglelo)
        vi(1, i) = loc_to_glo_nid(e)
      enddo

      nr  = nvi
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     $  vl, 0, vr, 0, key)

      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      if (wdsizi.eq.8) then
        do i = 1, nr
          call copyi4(eg(i), vi(2, i), 1)
        enddo
      else if (wdsizi.eq.4) then
        do i = 1, nr
          eg(i) = vi(2, i)
        enddo
      endif

      do i = 1, nr
         call find_lglel_ind(e, eg(i), nelt)
         call buf_to_bc_loc(e, cbl, bl, vi(2, i))
      enddo

 100  call err_chk(ierr, 'Error transferring .re2 boundary data$')
      return
      end
c-----------------------------------------------------------------------
      subroutine bc_to_buf(buf, cbl, bl, eg, f)

      include 'SIZE'
      include 'TOTAL'

      integer buf(30), eg, f
      character*3 cbl
      real bl(5)

      if (wdsizi.eq.8) then
        call icopy48(buf(1) , eg, 1) ! 1 - 2
        call icopy48(buf(3) ,  f, 1) ! 3 - 4
        call copy   (buf(5) , bl, 5) ! 5 -14
        call chcopy (buf(15),cbl, 3) !15
      else
        call icopy (buf(1), eg, 1) ! 1
        call icopy (buf(2),  f, 1) ! 2
        call copyX4(buf(3), bl, 5) ! 3 - 7
        call chcopy(buf(8),cbl, 3) ! 8
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_bc_big(cbl, bl, buf)

      include 'SIZE'
      include 'TOTAL'

      character*3 cbl(2*ldim, lelt)
      real bl(5, 2*ldim, lelt)

      integer e, f, eg, buf(30)

      if (wdsizi.eq.8) then
        call icopy84(eg, buf(1) ,1) ! 1 - 2
        call icopy84(f , buf(3), 1) ! 3 - 4
        call find_lglel_ind(e, eg, nelt)
        call copy  (bl(1, f, e), buf(5) , 5) ! 5 -14
        call chcopy(cbl(f, e)  , buf(15), 3) !15

        ! Integer assign connecting P element
        if (nelgt.ge.1000000.and.cbl(f, e).eq.'P  ')
     $   call copyi4(bl(1, f, e), buf(5), 1)

      else
        eg = buf(1)
        f  = buf(2)
        call find_lglel_ind(e, eg, nelt)
        call copy4r(bl(1, f, e), buf(3), 5)
        call chcopy(cbl(f, e)  , buf(8), 3)

        ! Integer assign of connecting periodic element
        if (nelgt.ge.1000000.and.cbl(f, e).eq.'P  ')
     $     bl(1, f, e) = buf(3)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine transfer_re2_curve(loc_to_glob_nid, lglelo, nelto)

      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = lelt*(2*ldim*(ldim - 1)))
      parameter(lrs = 15)
      parameter(li = 2*lrs + 1)

      ! character*1 ccurve(12, lelt)
      ! real curve(6, 12, lelt)
      integer loc_to_glob_nid(lelt), lglelo(lelt), nelto

      integer e, eg, f, cnt, nr, key, ierr, nf

      integer         vi(li, nrmax)
      common /ctmp1/  vi

      cnt = 0
      nedge = 2*ndim*(ndim - 1)
      do e = 1, nelto
        eg = lglelo(e)
        do f = 1, nedge
          cnt = cnt + 1
          vi(1, cnt) = loc_to_glob_nid(e)
          call curve_to_buf(vi(2, cnt), curve(1, f, e), ccurve(f, e),
     $      eg, f)
        enddo
      enddo

      nr = cnt
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     $  vl, 0, vr, 0, key)

      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      do i = 1,nr
         call buf_to_curve_big(curve, ccurve, vi(2, i))
      enddo

 100  call err_chk(ierr, 'Error transferring .re2 boundary data$')
      return
      end
c-----------------------------------------------------------------------
      subroutine transfer_re2_curve_v2(nvi, vi, loc_to_glo_nid, lglelo,
     $  nelto)
      include 'SIZE'
      include 'TOTAL'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      parameter(nrmax = 12*lelt)
      parameter(lrs = 8)
      parameter(li = 2*lrs + 1)

      integer nvi, vi(li, nrmax)
      integer loc_to_glo_nid(lelt), lglelo(lelt), nelto

      integer i, e, nr, key, ierr
      integer eg(nrmax)

      if (wdsizi.eq.8) then
        do i = 1, nvi
          call copyi4(eg(i), vi(2, i), 1)
        enddo
      else if (wdsizi.eq.4) then
        do i = 1, nvi
          eg(i) = vi(2, i)
        enddo
      endif

      do i = 1, nvi
        call find_sorted_ind(e, eg(i), nelto, lglelo)
        vi(1, i) = loc_to_glo_nid(e)
      enddo

      nr  = nvi
      key = 1
      call fgslib_crystal_tuple_transfer(cr_re2, nr, nrmax, vi, li,
     $  vl, 0, vr, 0, key)

      ierr = 0
      if (nr.gt.nrmax) then
        ierr = 1
        goto 100
      endif

      if (wdsizi.eq.8) then
        do i = 1, nr
          call copyi4(eg(i), vi(2, i), 1)
        enddo
      else if (wdsizi.eq.4) then
        do i = 1, nr
          eg(i) = vi(2, i)
        enddo
      endif

      do i = 1, nr
         call find_lglel_ind(e, eg(i), nelt)
         call buf_to_curve_loc(e, vi(2, i))
      enddo

 100  call err_chk(ierr, 'Error transferring .re2 curve data$')
      return
      end
c-----------------------------------------------------------------------
      subroutine curve_to_buf(buf, cv, ccv, eg, f)

      include 'SIZE'
      include 'TOTAL'

      integer buf(30), eg, f
      real cv(6)
      character*1 ccv

      if (wdsizi.eq.8) then
        call icopy48(buf(1) , eg, 1) ! 1 - 2
        call icopy48(buf(3) ,  f, 1) ! 3 - 4
        call copy   (buf(5) , cv, 5) ! 5 -14
        call chcopy (buf(15),ccv, 1) !15
      else
        call icopy (buf(1), eg, 1) ! 1
        call icopy (buf(2),  f, 1) ! 2
        call copyX4(buf(3), cv, 5) ! 3 - 7
        call chcopy(buf(8),ccv, 1) ! 8
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_curve_big(cv, ccv, buf)

      include 'SIZE'
      include 'TOTAL'

      character*1 ccv(12, lelt)
      real cv(6, 12, lelt)
      integer buf(30)

      integer e, f, eg

      if (wdsizi.eq.8) then
        call icopy84(eg, buf(1) ,1) ! 1 - 2
        call icopy84(f , buf(3), 1) ! 3 - 4
        call find_lglel_ind(e, eg, nelt)
        call copy  (cv(1, f, e), buf(5) , 5) ! 5 -14
        call chcopy(ccv(f, e)  , buf(15), 1) !15
      else
        eg = buf(1)
        f  = buf(2)
        call find_lglel_ind(e, eg, nelt)
        call copy4r(cv(1, f, e), buf(3), 5)
        call chcopy(ccv(f, e)  , buf(8), 1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine readat_big_v2
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'CTIMER'

      parameter(nrmax = 12*lelt) ! maximum number of records
      parameter(lrs   = 30)   ! record size: eg iside curve(5) ccurve
      parameter(li    = 2*lrs+2)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      common /ivrtx/ vertex((2**ldim),lelt)
      integer*8 vertex

      integer        vi(li, nrmax)
      common /ctmp1/ vi

      logical ifbswap
      integer loc_to_glo_nid(lelt),lglelo(lelt),nelto,nvi,ierr
      integer ibuf(2)
      integer iwork(lelt)
      real etimei, etimee, dur0, dur1

      ierr = 0

      ! Read the header to get nelgt, nelgv, ndim and then calculate
      ! nelt.
      call read_re2_hdr(ifbswap, .true.)

      ! Mesh is read by binning -- nelt or nelt + 1 elements in each
      ! rank (nelt + 1 on last nr ranks where nr is nelgt % np).
      nelt = nelgt / np
      do i = 1, mod(nelgt, np)
        if (np-i.eq.nid) nelt = nelt + 1
      enddo
      if (nelt.gt.lelt) then
        call exitti('nelt > lelt!$',nelt)
      endif

      call chkParam

#ifndef NOMPIIO
      call byte_open_mpi(re2fle, fh_re2, .true. ,ierr)
      call err_chk(ierr, ' Cannot open .re2 file!$')
#else
      call exitti('No serial support for big mesh read! P=$',np)
#endif

      if (nio.eq.0 .and. loglevel.gt.1) write(6,'(A)')
     $  ' partioning elements to MPI ranks'
      etimei = dnekclock_sync()

      dur0 = dnekclock_sync()
      call read_re2_mesh_v2(ifbswap)
      dur1 = dnekclock_sync() - dur0
      if (nio.eq.o .and. loglevel.gt.1) then
        write(6, *) 'done :: read coordinates: ', dur1
      endif

c     Distributed memory processor mapping
      if (np.gt.nelgt) then
         if(nio.eq.0) then
           write(6,1000) np,nelgt
 1000      format(2X,'ABORT: Too many processors (',I8
     $          ,') for to few elements (',I12,').'
     $          ,/,2X,'Aborting in readat_big_v2.')
         endif
         call exitt
      endif

      call get_vert_big_v2(vertex,loc_to_glo_nid)

      call fgslib_crystal_setup(cr_re2,nekcomm,np)

      ! transfer vertices
      nelto = nelt
      do i = 1, nelt
        lglelo(i) = lglel(i)
      enddo

      dur0 = dnekclock_sync()
      call transfer_vertices_v2(vertex,loc_to_glo_nid)
      dur1 = dnekclock_sync() - dur0
      if (nio.eq.0 .and. loglevel.gt.1) then
        write(6, *) 'done :: transfer vertices: ', dur1
      endif

      ! transfer coordinates
      dur0 = dnekclock_sync()
      call transfer_re2_mesh_v2(loc_to_glo_nid, lglelo, nelto)
      dur1 = dnekclock_sync() - dur0
      if (nio.eq.0 .and. loglevel.gt.1) then
        write(6, *) 'done :: transfer coordinates: ', dur1
      endif

      ! read and transfer curve sides
      dur0 = dnekclock_sync()
      call read_re2_curve_v2(nvi, vi, ifbswap)
      dur1 = dnekclock_sync() - dur0
      if (nio.eq.0 .and. loglevel.gt.1) then
        write(6, *) 'done :: read curved sides: ', dur1
      endif

      dur0 = dnekclock_sync()
      call transfer_re2_curve_v2(nvi, vi, loc_to_glo_nid, lglelo,
     $  nelto)
      dur1 = dnekclock_sync() - dur0
      if (nio.eq.0 .and. loglevel.gt.1) then
        write(6, *) 'done :: transfer curved sides: ', dur1
      endif

      ! transfer bcs
      dur0 = dnekclock_sync()

                  ibc = 2
      if (ifflow) ibc = 1

                  nfldt = 1
      if (ifheat) nfldt = 2+npscal
      if (ifmhd ) nfldt = 2+npscal+1

      ! first field to read
      if (param(33).gt.0) ibc = int(param(33))

      ! number of fields to read
      if (param(32).gt.0) then
        nfldt = ibc + int(param(32)) - 1
        nfldt = max(nfldt,1)
        if (nelgt.gt.nelgv) nfldt = max(nfldt,2)
      endif

      call blank(cbc,3*size(cbc))
      call rzero(bc ,size(bc))

      do ifield = ibc, nfldt
        call read_re2_bc_v2(nvi, vi, cbc(1,1,ifield), bc(1,1,1,ifield),
     $    ifbswap)

        call transfer_re2_bc_v2(nvi, vi, cbc(1,1,ifield),
     $    bc(1,1,1,ifield), loc_to_glo_nid, lglelo, nelto)
      enddo

      dur1 = dnekclock_sync() - dur0
      if (nio.eq.0 .and. loglevel.gt.1) then
        write(6, *) 'done :: read and transfer bcs: ', dur1
      endif

      call fgslib_crystal_free(cr_re2)

#ifndef NOMPIIO
      call byte_close_mpi(fh_re2,ierr)
#endif

#ifdef DPROCMAP
      call dProcmapInit()
      do i = 1,nelt
         ieg = lglel(i)
         if (ieg.lt.1 .or. ieg.gt.nelgt)
     $      call exitti('invalid ieg!$',ieg)
         ibuf(1) = i
         ibuf(2) = nid
         call dProcmapPut(ibuf,2,0,ieg)
      enddo
#else
      call izero(gllel,nelgt)
      do i = 1,nelt
         ieg = lglel(i)
         if (ieg.lt.1 .or. ieg.gt.nelgt)
     $      call exitti('invalid ieg!$',ieg)
         gllnid(ieg) = nid
         gllel(ieg) = i
      enddo
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllnid(k),iwork,'+  ',m)
         if (m.gt.0) call igop(gllel(k) ,iwork,'+  ',m)
         k = k+m
      enddo
#endif

      etimee = dnekclock_sync() - etimei
      if (nio.eq.0 .and. loglevel.gt.1) then
        write(6, *) 'done :: partitioning: ', etimee
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine read_re2_mesh_v2(ifbswap)
      include 'SIZE'
      include 'TOTAL'

      parameter(nrmax = lelt)             ! maximum number of records
      parameter(lrs   = 1+ldim*(2**ldim)) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      integer         bufr(li-2, nrmax)
      common /scrns/  bufr

      logical ifbswap
      integer e, eg, ind(nrmax), ierr

      integer*8 lre2off_b, dtmp8, nrg

      nrg       = nelgt
      nr        = nelt
      irankoff  = igl_running_sum(nr) - nr
      dtmp8     = irankoff
      re2off_b  = 84 ! set initial offset (hdr + endian)
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      ! read coordinates from file
      ierr = 0
      nwds4r = nr*lrs4
      call byte_set_view(lre2off_b,fh_re2)
      call byte_read_mpi(bufr,nwds4r,-1,fh_re2,ierr)
      re2off_b = re2off_b + nrg*4*lrs4
      if (ierr.gt.0) goto 100

      if (nio.eq.0) write(6,*) 'reading mesh '

      ! pack buffer
      do e = 1, nr
         lglel(e)= irankoff + e ! elements are stored in global order
         jj      = (e-1)*lrs4 + 1
         call buf_to_xyz(bufr(jj,1),e,ifbswap,ierr)
      enddo

 100  call err_chk(ierr,'Error reading .re2 mesh$')

      return
      end
c-----------------------------------------------------------------------
      subroutine get_dest_proc(p, eg, np, nelgt)
        integer p, eg, np, nelgt
        integer nlt, nr, nstar

        nlt = nelgt / np
        nr = nelgt - nlt * np
        nstar = nlt * (np - nr)

        if (eg.le.nstar) then
          p = (eg - 1) / nlt
        else
          p = (eg - nstar - 1) / (nlt + 1) + np - nr
        endif
        return
      end
c-----------------------------------------------------------------------
